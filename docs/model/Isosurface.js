import * as THREE from '../external/three/three.module.js';
// import { MarchCubes } from '../compiled/MarchCubes.js';

import { groups } from '../state/store.js';

import { MarchingCubesWrapper, MarchingCubesBackend } from './MarchingCubesWrapper.js';
import { applyTransparency } from '../utils/TransparencyPolicy.js';
import { makeFractionalBoundsClippingPlanes, ColormapLut } from './Plane.js';


// Transparency flags (transparent, depthWrite, renderOrder) are owned by the
// rendering pipeline policy ('isosurface' in render/pipeline/ForwardPipeline.js)
// and applied via applyIsosurfaceMaterialSettings below.
export let surface_options = {
    //transmission: 0.95,
    side: THREE.DoubleSide,
    opacity:0.6,
  }



export let defaultPosColor = new THREE.Color(0x33aaff);
const _white = new THREE.Color(0xffffff);
export let defaultNegColor = new THREE.Color(0xff3333);
export let isosurfaceTriangleSortingEnabled = true;


const _sortCameraPosition = new THREE.Vector3();
const _sortCameraQuaternion = new THREE.Quaternion();
const _lastSortCameraPosition = new THREE.Vector3();
const _lastSortCameraQuaternion = new THREE.Quaternion();
let _hasLastSortCameraState = false;

function clampOpacity(opacity) {
    if (!Number.isFinite(opacity)) return surface_options.opacity;
    return Math.max(0, Math.min(1, opacity));
}

export function getIsosurfaceMaterialSettings() {
    return {
        positiveColor: `#${defaultPosColor.getHexString()}`,
        negativeColor: `#${defaultNegColor.getHexString()}`,
        opacity: surface_options.opacity
    };
}

export function setIsosurfaceMaterialSettings(settings = {}) {
    if (settings.positiveColor !== undefined) {
        defaultPosColor.set(settings.positiveColor);
    }
    if (settings.negativeColor !== undefined) {
        defaultNegColor.set(settings.negativeColor);
    }
    if (settings.opacity !== undefined) {
        surface_options.opacity = clampOpacity(settings.opacity);
    }
}

// ---------------------------------------------------------------------------
//  Isosurface smoothing — one active method, GLOBAL like the material settings
//  above, with each method's parameters remembered separately so switching
//  back restores them. The marching-cubes wrapper (setSmoothing in
//  model/MarchingCubesWrapper.js) consumes these ids and param keys verbatim.
//
//  'mesh' methods post-process the triangle mesh ('project' snaps vertices
//  back onto the exact isovalue afterwards); 'field' methods alter the grid
//  before meshing. The ray/path tracers render the raw field, so mesh methods
//  affect only the rasterised surface.
// ---------------------------------------------------------------------------

/**
 * @typedef {object} SmoothingParamSpec
 * @property {string} key
 * @property {string} label
 * @property {'int'|'float'|'bool'|'choice'} type
 * @property {number} [min]
 * @property {number} [max]
 * @property {number} [step]
 * @property {number|boolean|string} default
 * @property {ReadonlyArray<{value: string, label: string}>} [choices]
 */

/**
 * @typedef {object} SmoothingMethodSpec
 * @property {string} id
 * @property {string} label
 * @property {'mesh'|'field'|null} domain
 * @property {ReadonlyArray<SmoothingParamSpec>} params
 */

/**
 * @typedef {object} IsosurfaceSmoothingSettings
 * @property {string} method active method id ('off' for none)
 * @property {Record<string, Record<string, number|boolean|string>>} params per-method values
 */

const projectParam = { key: 'project', label: 'Project onto isovalue', type: 'bool', default: true };

/** @type {ReadonlyArray<SmoothingMethodSpec>} */
export const SMOOTHING_METHODS = Object.freeze(/** @type {SmoothingMethodSpec[]} */ ([
    { id: 'off', label: 'Off', domain: null, params: [] },
    {
        id: 'laplacian', label: 'Laplacian', domain: 'mesh',
        params: [
            { key: 'iterations', label: 'Iterations', type: 'int', min: 1, max: 50, step: 1, default: 5 },
            { key: 'lambda', label: 'λ', type: 'float', min: 0.05, max: 1, step: 0.05, default: 0.5 },
            projectParam,
        ],
    },
    {
        id: 'taubin', label: 'Taubin λ|μ', domain: 'mesh',
        params: [
            { key: 'iterations', label: 'Iterations', type: 'int', min: 1, max: 100, step: 1, default: 10 },
            { key: 'lambda', label: 'λ', type: 'float', min: 0.05, max: 0.9, step: 0.05, default: 0.5 },
            { key: 'passband', label: 'Pass-band', type: 'float', min: 0.01, max: 0.3, step: 0.01, default: 0.1 },
            projectParam,
        ],
    },
    {
        id: 'hc', label: 'HC-Laplacian', domain: 'mesh',
        params: [
            { key: 'iterations', label: 'Iterations', type: 'int', min: 1, max: 50, step: 1, default: 10 },
            { key: 'alpha', label: 'α', type: 'float', min: 0, max: 1, step: 0.05, default: 0.1 },
            { key: 'beta', label: 'β', type: 'float', min: 0, max: 1, step: 0.05, default: 0.5 },
            projectParam,
        ],
    },
    {
        id: 'loop', label: 'Loop subdivision', domain: 'mesh',
        params: [
            { key: 'levels', label: 'Levels', type: 'int', min: 1, max: 3, step: 1, default: 1 },
            projectParam,
        ],
    },
    {
        id: 'catmullClark', label: 'Catmull–Clark', domain: 'mesh',
        params: [
            { key: 'levels', label: 'Levels', type: 'int', min: 1, max: 2, step: 1, default: 1 },
            projectParam,
        ],
    },
    {
        id: 'tricubic', label: 'Tricubic refinement', domain: 'field',
        params: [
            { key: 'factor', label: 'Factor', type: 'int', min: 2, max: 4, step: 1, default: 2 },
            {
                key: 'kernel', label: 'Kernel', type: 'choice', default: 'catmullRom',
                choices: [
                    { value: 'catmullRom', label: 'Catmull-Rom' },
                    { value: 'bspline', label: 'B-spline' },
                ],
            },
        ],
    },
    {
        id: 'gaussian', label: 'Gaussian pre-filter', domain: 'field',
        params: [
            { key: 'sigma', label: 'σ (voxels)', type: 'float', min: 0.3, max: 3, step: 0.1, default: 1.0 },
        ],
    },
]).map((method) => Object.freeze({
    ...method,
    params: Object.freeze(method.params.map((p) => Object.freeze({
        ...p, ...(p.choices ? { choices: Object.freeze(p.choices.map((c) => Object.freeze({ ...c }))) } : {}),
    }))),
})));

/** @param {string} id @returns {SmoothingMethodSpec | undefined} */
function findSmoothingMethod(id) {
    return SMOOTHING_METHODS.find((m) => m.id === id);
}

/** Coerce one value to its spec (clamped, stepped for ints); null if unusable. */
function sanitizeSmoothingParam(spec, value) {
    switch (spec.type) {
        case 'bool':
            return typeof value === 'boolean' ? value : null;
        case 'choice':
            return spec.choices?.some((c) => c.value === value) ? value : null;
        default: {
            let v = Number(value);
            if (typeof value === 'boolean' || value === null || value === '' || !Number.isFinite(v)) return null;
            if (spec.type === 'int') v = Math.round(v);
            return Math.max(spec.min ?? -Infinity, Math.min(spec.max ?? Infinity, v));
        }
    }
}

const _smoothing = {
    method: 'off',
    /** @type {Record<string, Record<string, number|boolean|string>>} */
    params: Object.fromEntries(SMOOTHING_METHODS.map((m) =>
        [m.id, Object.fromEntries(m.params.map((p) => [p.key, p.default]))])),
};

/** @returns {IsosurfaceSmoothingSettings} a deep copy of the live settings */
export function getIsosurfaceSmoothingSettings() {
    return {
        method: _smoothing.method,
        params: Object.fromEntries(Object.entries(_smoothing.params).map(([id, p]) => [id, { ...p }])),
    };
}

/**
 * Merge a (partial) setting: an unknown method becomes 'off', per-method
 * params are clamped to their spec, unknown methods/keys are ignored.
 * @param {{method?: string, params?: Record<string, Record<string, any>>}} [settings]
 */
export function setIsosurfaceSmoothingSettings(settings = {}) {
    if (settings.method !== undefined) {
        _smoothing.method = findSmoothingMethod(settings.method) ? settings.method : 'off';
    }
    const params = settings.params;
    if (!params || typeof params !== 'object') return;
    for (const [id, values] of Object.entries(params)) {
        const method = findSmoothingMethod(id);
        if (!method || !values || typeof values !== 'object') continue;
        for (const spec of method.params) {
            if (!(spec.key in values)) continue;
            const v = sanitizeSmoothingParam(spec, values[spec.key]);
            if (v !== null) _smoothing.params[id][spec.key] = v;
        }
    }
}

export function setIsosurfaceTriangleSortingEnabled(enabled) {
    isosurfaceTriangleSortingEnabled = Boolean(enabled);
}

export function getIsosurfaceTriangleSortingEnabled() {
    // return isosurfaceTriangleSortingEnabled;
    // Disabled as the rendering engine fixes the artifacts this was meant to
    // solve: below alpha 1 the order-independent pipelines (depth peeling /
    // WBOIT) blend the surface correctly, and AT alpha 1 the transparency
    // policy now marks it genuinely opaque (depth writes on, see
    // render/pipeline/ForwardPipeline.js 'isosurface'), so the depth buffer —
    // not the marching-cubes vertex order — decides what is in front.
    return false;
}

// ---------------------------------------------------------------------------
//  Colour-by: an isosurface of field A coloured per vertex by field B
//
//  Field.colorBy = { field: Field, colormap: string, min: number, max: number }
//  (or null). The marching-cubes backend interpolates B onto every vertex
//  with the same edge parameter as the vertex position
//  (model/MarchingCubesWrapper.js setColorField), the raw values are kept on
//  the geometry (userData.colorByValues) so a colormap/range edit only
//  remaps, and the mapped RGB (userData.colorByRGB) goes into the RGBA vertex
//  `color` attribute that render/FocusRegionModule.js shares for its alpha
//  (userData.focusAlpha). composeIsosurfaceVertexColors() is the single
//  writer of that attribute, so whichever of the two runs first, the result
//  is the same.
// ---------------------------------------------------------------------------

/**
 * The field's colour-by setting if it can be honoured: a colour field with
 * values on the same grid and a finite range. Anything else reads as none.
 *
 * @param {any} field a Field (colorBy may not be declared on older builds)
 * @returns {{field: any, colormap: string, min: number, max: number} | null}
 */
export function resolveColorBy(field) {
    const colorBy = field?.colorBy;
    const source = colorBy?.field;
    if (!source?.values || !field) return null;
    if (source.nx !== field.nx || source.ny !== field.ny || source.nz !== field.nz) return null;
    const min = Number(colorBy.min);
    const max = Number(colorBy.max);
    if (!Number.isFinite(min) || !Number.isFinite(max)) return null;
    return { field: source, colormap: colorBy.colormap || 'bgyor', min, max };
}

/**
 * Whether `candidate` can colour `field`'s isosurface: loaded and on the same
 * grid. (A field may colour itself — that is just a gradient along the
 * isovalue, harmless and occasionally a useful check.)
 */
export function canColorFieldBy(field, candidate) {
    return Boolean(field && candidate?.values)
        && candidate.nx === field.nx && candidate.ny === field.ny && candidate.nz === field.nz;
}

/**
 * Map per-vertex values through a colormap into linear RGB triplets.
 * @param {ArrayLike<number>} values
 * @param {{colormap: string, min: number, max: number}} colorBy
 * @returns {Float32Array}
 */
export function mapValuesToColors(values, colorBy) {
    const lut = new ColormapLut(colorBy.colormap).setMin(colorBy.min).setMax(colorBy.max);
    const out = new Float32Array(values.length * 3);
    for (let i = 0; i < values.length; i++) {
        const v = values[i];
        // NaN (never produced by marching cubes, but a colour field may hold
        // them) maps to the low end rather than poisoning the attribute.
        const c = lut.getColor(Number.isFinite(v) ? v : colorBy.min);
        out[i * 3] = c.r;
        out[i * 3 + 1] = c.g;
        out[i * 3 + 2] = c.b;
    }
    return out;
}

/**
 * Write (or remove) an isosurface mesh's RGBA vertex `color` attribute from
 * its two independent inputs on the geometry:
 *   userData.colorByRGB  Float32Array(count*3) | null — colour-by RGB (else 1)
 *   userData.focusAlpha  Float32Array(count)   | null — focus fade (else 1)
 * With neither, the attribute is removed and vertex colours switched off.
 * Transparency is NOT touched here: the focus code owns the
 * needsTransparency decision (its alpha), colour-by is opaque RGB.
 *
 * @param {any} mesh
 */
export function composeIsosurfaceVertexColors(mesh) {
    const geometry = mesh?.geometry;
    const material = mesh?.material;
    if (!geometry || !material) return;
    const position = geometry.getAttribute('position');
    const count = position?.count ?? 0;
    const rgb = geometry.userData?.colorByRGB;
    const alpha = geometry.userData?.focusAlpha;
    const hasRgb = Boolean(rgb) && rgb.length === count * 3 && count > 0;
    const hasAlpha = Boolean(alpha) && alpha.length === count && count > 0;
    if (!hasRgb && !hasAlpha) {
        if (geometry.getAttribute('color')) geometry.deleteAttribute('color');
        if (material.vertexColors) {
            material.vertexColors = false;
            material.needsUpdate = true;
        }
        return;
    }
    let color = geometry.getAttribute('color');
    if (!color || color.itemSize !== 4 || color.count !== count) {
        color = new THREE.BufferAttribute(new Float32Array(count * 4), 4);
        geometry.setAttribute('color', color);
    }
    const dst = /** @type {Float32Array} */ (color.array);
    for (let i = 0; i < count; i++) {
        dst[i * 4] = hasRgb ? rgb[i * 3] : 1;
        dst[i * 4 + 1] = hasRgb ? rgb[i * 3 + 1] : 1;
        dst[i * 4 + 2] = hasRgb ? rgb[i * 3 + 2] : 1;
        dst[i * 4 + 3] = hasAlpha ? alpha[i] : 1;
    }
    color.needsUpdate = true;
    if (!material.vertexColors) {
        material.vertexColors = true;
        material.needsUpdate = true;
    }
}

/**
 * Re-colour the live isosurface after its field's colorBy changed.
 *
 * A colormap or min/max edit only remaps the per-vertex values already on
 * the geometry. A different colour FIELD needs new values, which only a
 * marching-cubes pass produces, so then this returns false and the caller
 * rebuilds (render/index.js updateField). Callers still requestRender().
 *
 * @param {any} [isosurface] defaults to groups.isosurfaceGroup
 * @returns {boolean} true if the colours are now current without a rebuild
 */
export function refreshIsosurfaceColors(isosurface = groups.isosurfaceGroup) {
    if (!isosurface?.applyColorBy) return false;
    return isosurface.applyColorBy();
}

export function applyIsosurfaceMaterialSettings(isosurface, settings = {}) {
    if (!isosurface || !isosurface.meshes) return;

    const { positiveColor, negativeColor, opacity } = settings;

    const applyToMesh = (mesh, color) => {
        const material = mesh?.material;
        if (!material) return;
        if (color !== undefined) {
            // A colour-by surface takes its colour from the vertex attribute;
            // a tinted material would multiply into it. The lobe colour is
            // kept in defaultPos/NegColor and restored when colour-by ends.
            material.color.set(mesh.geometry?.userData?.colorByRGB ? 0xffffff : color);
        }
        if (opacity !== undefined) {
            material.opacity = clampOpacity(opacity);
        }
        // Keep a focus-region vertex fade (render/FocusRegionModule.js
        // applyFocusToField) in the blended pass across an opacity edit.
        applyTransparency(material, {
            kind: 'isosurface', opacity: material.opacity,
            needsTransparency: !!material.userData?.transparencySpec?.needsTransparency, mesh,
        });
    };

    applyToMesh(isosurface.meshes.positive, positiveColor);
    applyToMesh(isosurface.meshes.negative, negativeColor);
    // The periodic boundary copies share these materials but carry their own
    // renderOrder/visibility, which the transparency policy just rewrote.
    isosurface._syncImageState?.();
}

export function applyMaterialSettingsToStoredIsosurfaces(isosurfaceGroup, settings = {}) {
    if (!isosurfaceGroup) return;

    const applyOne = (entry) => {
        if (!entry) return;
        if (entry.meshes?.positive || entry.meshes?.negative) {
            applyIsosurfaceMaterialSettings(entry, settings);
            return;
        }
        if (entry.traverse) {
            entry.traverse((child) => {
                if (child?.meshes?.positive || child?.meshes?.negative) {
                    applyIsosurfaceMaterialSettings(child, settings);
                }
            });
        }
    };

    if (Array.isArray(isosurfaceGroup)) {
        isosurfaceGroup.forEach(applyOne);
    } else if (isosurfaceGroup instanceof Set) {
        isosurfaceGroup.forEach(applyOne);
    } else if (isosurfaceGroup instanceof Map) {
        isosurfaceGroup.forEach((value) => applyOne(value));
    } else {
        applyOne(isosurfaceGroup);
    }
}

export function updateStoredIsosurfaceRenderOrder(camera, isosurfaceGroup) {
    if (!camera || !isosurfaceGroup) return;

    camera.updateMatrixWorld(true);
    _sortCameraPosition.setFromMatrixPosition(camera.matrixWorld);
    camera.getWorldQuaternion(_sortCameraQuaternion);

    if (
        _hasLastSortCameraState
        && _lastSortCameraPosition.equals(_sortCameraPosition)
        && _lastSortCameraQuaternion.equals(_sortCameraQuaternion)
    ) {
        return;
    }

    _hasLastSortCameraState = true;
    _lastSortCameraPosition.copy(_sortCameraPosition);
    _lastSortCameraQuaternion.copy(_sortCameraQuaternion);

    const applyOne = (entry) => {
        if (!entry) return;
        if (typeof entry.sortTrianglesByCameraDistance === 'function') {
            entry.sortTrianglesByCameraDistance(_sortCameraPosition);
            return;
        }
        if (entry.traverse) {
            entry.traverse((child) => {
                if (typeof child?.sortTrianglesByCameraDistance === 'function') {
                    child.sortTrianglesByCameraDistance(_sortCameraPosition);
                }
            });
        }
    };

    if (Array.isArray(isosurfaceGroup)) {
        isosurfaceGroup.forEach(applyOne);
    } else if (isosurfaceGroup instanceof Set) {
        isosurfaceGroup.forEach(applyOne);
    } else if (isosurfaceGroup instanceof Map) {
        isosurfaceGroup.forEach((value) => applyOne(value));
    } else {
        applyOne(isosurfaceGroup);
    }
}



// ---------------------------------------------------------------------------
//  Periodic display boundary (VESTA-style "Active Cell Boundary")
//
//  general.periodicBounds gives a per-axis fractional [min, max] display
//  region, and render/LatticeModule.js draws every periodic image of every
//  atom that lands inside it. A volumetric field is periodic in exactly the
//  same way, so it follows the boundary the same way — it is just expressed
//  differently: instead of emitting extra atoms, the field is DRAWN AGAIN in
//  every cell the region reaches (the marching-cubes mesh is translated by
//  whole lattice vectors — same geometry, same material, one extra draw call)
//  and every copy is CLIPPED to the region, so a boundary that stops
//  part-way through a cell cuts the surface there instead of showing a whole
//  extra cell of it.
// ---------------------------------------------------------------------------

/** @type {[number, number][]} */
const UNIT_BOUNDS = [[0, 1], [0, 1], [0, 1]];
// Bounds are user-typed, so a value a hair over an integer (1.0000001) must
// not conjure a whole extra cell of field.
const BOUND_EPS = 1e-6;
// Safety net for a restored/shared state with wild bounds: the panel itself
// clamps to +/-2 cells (5 per axis), and 5^3 copies of one mesh is already a
// lot of geometry to push per frame.
const MAX_IMAGES_PER_AXIS = 5;

/** True for the plain unit cell [0,1] on every axis — the default, in which
 *  the field is exactly one copy and needs no clipping at all. */
function isUnitBounds(bounds) {
    return bounds.every(([lo, hi]) => Math.abs(lo) < BOUND_EPS && Math.abs(hi - 1) < BOUND_EPS);
}

/** Integer cell translations n whose own cell [n, n+1] overlaps [lo, hi].
 *  [0,1] -> [0] (today's single copy); [0,1.2] -> [0,1]; [-0.5,1] -> [-1,0]. */
function axisImageRange([lo, hi]) {
    const first = Math.floor(lo + BOUND_EPS);
    // max(): a zero-thickness region (lo === hi) still resolves to one cell,
    // which the clipping then reduces to nothing — better than no mesh at all.
    const last = Math.min(Math.max(first, Math.ceil(hi - BOUND_EPS) - 1), first + MAX_IMAGES_PER_AXIS - 1);
    const out = [];
    for (let n = first; n <= last; n++) out.push(n);
    return out;
}

/** Every integer cell translation [i,j,k] the display boundary reaches. */
function boundsImageOffsets(bounds) {
    const [ri, rj, rk] = bounds.map(axisImageRange);
    const out = [];
    for (const i of ri) for (const j of rj) for (const k of rk) out.push([i, j, k]);
    return out;
}

export class Isosurface extends THREE.Group{

    constructor(field) {
        super();
        this.field = field;
        /** @type {string} */
        this.backend = MarchingCubesBackend.WASM;
        this.lastCameraPosition = new THREE.Vector3();

        this.marchingCubes = new MarchingCubesWrapper(field, this.backend);
        /** The colour-by field the backend currently interpolates (null = none).
         *  @type {any} */
        this._colorSourceField = null;

        /** Extra meshes drawing this field in the other cells the periodic
         *  display boundary reaches. They SHARE the positive/negative meshes'
         *  geometry and material — only the cell translation differs. */
        this._imageMeshes = [];
        /** @type {[number, number][]} the boundary these copies were built for */
        this._periodicBounds = UNIT_BOUNDS;
        /** The structure lattice the boundary is expressed in (rows = vectors,
         *  world units). Null = this field's own grid cell. The two differ once
         *  a supercell is built after the field was loaded: the field is NOT
         *  duplicated, it keeps its original sub-cell, and the copies and the
         *  clipping box have to repeat that sub-cell (gaps included) with the
         *  structure's lattice, not tile the grid cell contiguously. */
        this._lattice = null;

        this.addMeshes();

        this.matrixAutoUpdate = false;
        const transform_cell = new THREE.Matrix4();
        transform_cell.set(
            this.field.voxel[0][0], this.field.voxel[1][0], this.field.voxel[2][0], 0,
            this.field.voxel[0][1], this.field.voxel[1][1], this.field.voxel[2][1], 0,
            this.field.voxel[0][2], this.field.voxel[1][2], this.field.voxel[2][2], 0,
            0, 0, 0, 1
        );
        transform_cell.scale(new THREE.Vector3(this.field.nx, this.field.ny, this.field.nz));
        this.applyMatrix4(transform_cell);
    }

    sortTrianglesByCameraDistance(cameraPosition) {
        if (!cameraPosition || !this.marchingCubes) return;

        // Transform camera position into the isosurface's local coordinate system
        this.updateMatrixWorld(true);
        this.lastCameraPosition.copy(cameraPosition);
        const localCameraPosition = cameraPosition.clone();
        const inverseMatrix = new THREE.Matrix4().copy(this.matrixWorld).invert();
        localCameraPosition.applyMatrix4(inverseMatrix);

        for (const meshKey of ['positive', 'negative']) {
            const mesh = this.meshes?.[meshKey];
            if (!mesh?.geometry) continue;

            const positionAttr = mesh.geometry.getAttribute('position');
            if (!positionAttr?.array || positionAttr.count < 3) continue;
            const normalAttr = mesh.geometry.getAttribute('normal');
            const colorAttr = mesh.geometry.getAttribute('color');
            const userData = mesh.geometry.userData ?? {};

            // Every per-vertex array rides the same permutation, or the
            // colour-by / focus colours would detach from their vertices.
            this.marchingCubes.sortVerticesToCamera(localCameraPosition, positionAttr.array, normalAttr?.array,
                colorAttr?.array, userData.colorByValues, userData.colorByRGB, userData.focusAlpha);

            positionAttr.needsUpdate = true;
            if (normalAttr) {
                normalAttr.needsUpdate = true;
            }
            if (colorAttr) colorAttr.needsUpdate = true;
        }
    }

    addMeshes() {
        let material_options = {};
        Object.assign(material_options, surface_options);
        material_options.color = defaultPosColor;
        const materialPos = new THREE.MeshPhysicalMaterial(material_options);
        material_options.color = defaultNegColor;
        const materialNeg = new THREE.MeshPhysicalMaterial(material_options);
        const positiveGeom = new THREE.BufferGeometry();
        const negativeGeom = new THREE.BufferGeometry();
        
        /** @type {{positive: any, negative: any}} */
        this.meshes = {
            positive: new THREE.Mesh(positiveGeom, materialPos),
            negative: new THREE.Mesh(negativeGeom, materialNeg)
        };
        this.meshes.positive.name = 'isosurface_pos';
        this.meshes.negative.name = 'isosurface_neg';

        // Applies the pipeline transparency policy too (incl. renderOrder).
        this.applyMaterialSettings(getIsosurfaceMaterialSettings());

        this.add(this.meshes.positive);
        this.add(this.meshes.negative);
    }

    applyMaterialSettings(settings = {}) {
        applyIsosurfaceMaterialSettings(this, settings);
    }

    get positiveMesh() {
        return this.meshes.positive;
    }

    get negativeMesh() {
        return this.meshes.negative;
    }

    set isovalue(value) {
        this.field.isovalue = value;
    }

    _replaceGeometry(meshKey, vertices, normals, values = null) {
        const tmpGeom = new THREE.BufferGeometry();
        tmpGeom.setAttribute('position', new THREE.BufferAttribute(vertices, 3));
        tmpGeom.setAttribute('normal', new THREE.BufferAttribute(normals, 3));
        // Colour-field values per vertex (colour-by), and which field they
        // came from, so applyColorBy() can tell a remap from a stale set.
        tmpGeom.userData.colorByValues = values ?? null;
        tmpGeom.userData.colorBySource = values ? this._colorSourceField : null;

        const merged = tmpGeom; // mergeVertices(tmpGeom); // merging vertices might be a good idea, but currently way more expensive
        if (this.backend === MarchingCubesBackend.THREE) {
            // the JS/THREE gets the normal wrong somehow, so we recompute it here.
            // The WASM backend computes correct normals, so we can skip this step for it.
            merged.computeVertexNormals();
        }
        //merged.computeBoundingSphere();

        this.meshes[meshKey].geometry = merged;
        this._syncImageState(); // the copies draw this same geometry
    }

    refreshGeometry(field_key, vertices, normals, values = null) {
        if (field_key == "positive") {
            if (vertices === undefined || normals === undefined) {
                const data = this.marchingCubes.getVertices();
                vertices = data.vertices;
                normals = data.normals;
                values = data.values;
            }
            this._replaceGeometry('positive', vertices, normals, values);
            this.meshes.positive.needsUpdate = true;
            this.add(this.meshes.positive);
        }
        else if (field_key == "negative") {
            if (vertices === undefined || normals === undefined) {
                const data = this.marchingCubes.getVertices();
                vertices = data.vertices;
                normals = data.normals;
                values = data.values;
            }
            this._replaceGeometry('negative', vertices, normals, values);
            this.meshes.negative.needsUpdate = true;
            this.add(this.meshes.negative);
        }
    }

    /**
     * A representative first isosurface level for this field, computed by the
     * marching-cubes backend from the values it already holds.
     *
     * @param {number} [fraction] target fraction of the cell inside the surface
     * @returns {number | null} null if the backend cannot supply one
     */
    defaultIsoValue(fraction = 0.05) {
        return this.marchingCubes ? this.marchingCubes.defaultIsoValue(fraction) : null;
    }

    /**
     * Follow the periodic display boundary: draw this field in every cell the
     * boundary reaches and clip every copy to it.
     *
     * Idempotent and cheap — no marching cubes rerun, the copies share the
     * base meshes' geometry and material — so it is safe to call whenever the
     * boundary, the field or the pipeline changes.
     *
     * @param {[number, number][]} [bounds] per-axis [min, max] in fractional
     *   coordinates OF THE STRUCTURE LATTICE, as render/LatticeModule.js
     *   normalizePeriodicBounds() returns. Defaults to the plain unit cell.
     * @param {number[][]|null} [lattice] the structure lattice those fractions
     *   refer to (rows = vectors, world units). Omitted: this field's own cell.
     */
    setPeriodicBounds(bounds = UNIT_BOUNDS, lattice = null) {
        const safe = /** @type {[number, number][]} */ (
            Array.isArray(bounds) && bounds.length === 3
                && bounds.every((b) => Array.isArray(b) && b.length === 2 && b.every(Number.isFinite))
                ? bounds.map(([lo, hi]) => (lo <= hi ? [lo, hi] : [hi, lo]))
                : UNIT_BOUNDS);
        const validLattice = Array.isArray(lattice) && lattice.length === 3
            && lattice.every((row) => Array.isArray(row) && row.length === 3 && row.every(Number.isFinite));
        this._periodicBounds = safe;
        this._lattice = validLattice ? lattice.map((row) => [...row]) : null;
        this._syncImageMeshes(boundsImageOffsets(safe));
        this._applyBoundsClipping(safe);
    }

    /** The lattice the boundary repeats the field with: the structure's when
     *  known, else this field's own grid cell. */
    _boundaryVectors() {
        return this._lattice ?? this._cellVectors();
    }

    /** True when the boundary lattice is this field's own grid cell (within
     *  rounding), i.e. the field fills exactly one boundary cell. */
    _latticeIsOwnCell() {
        if (!this._lattice) return true;
        const own = this._cellVectors();
        return own.every((vec, i) => vec.every((v, k) => Math.abs(v - this._lattice[i][k]) < 1e-6));
    }

    /** A whole-lattice translation [i,j,k] as a position in this group's local
     *  space. Local space is the grid's fractional coordinates, so the
     *  translation is the integer offset itself when the boundary lattice IS
     *  the grid cell, and the (linear) inverse of the group matrix applied to
     *  i*a + j*b + k*c otherwise. */
    _imagePosition([i, j, k]) {
        if (this._latticeIsOwnCell()) return new THREE.Vector3(i, j, k);
        const [a, b, c] = this._lattice;
        const world = new THREE.Vector3(
            i * a[0] + j * b[0] + k * c[0],
            i * a[1] + j * b[1] + k * c[1],
            i * a[2] + j * b[2] + k * c[2]);
        const toLocal = new THREE.Matrix3().setFromMatrix4(new THREE.Matrix4().copy(this.matrix).invert());
        return world.applyMatrix3(toLocal);
    }

    /** The cell this field is drawn in, read back out of the group's own
     *  voxel*dims matrix (its columns ARE the lattice vectors), so the copies
     *  and the clipping box use exactly the basis the surface is drawn in. */
    _cellVectors() {
        const e = this.matrix.elements;
        return [
            [e[0], e[1], e[2]],
            [e[4], e[5], e[6]],
            [e[8], e[9], e[10]],
        ];
    }

    /** One mesh per (cell translation x lobe). The FIRST translation is given
     *  to the base meshes themselves — the group's local space is fractional
     *  coordinates, so a whole-cell translation is just an integer position —
     *  and the rest become shared-geometry copies. */
    _syncImageMeshes(offsets) {
        const [first, ...rest] = offsets.length ? offsets : [[0, 0, 0]];
        const firstPosition = this._imagePosition(first);
        for (const key of ['positive', 'negative']) {
            this.meshes?.[key]?.position.copy(firstPosition);
        }

        // Rebuild the copy list only when the cells themselves (or the lattice
        // they are translated by) changed; a slider drag inside one cell then
        // costs nothing but the clip update.
        const latticeKey = (this._lattice ?? []).flat().map((v) => v.toFixed(6)).join(',');
        const key = `${rest.map((o) => o.join(',')).join(';')}|${latticeKey}`;
        if (key !== this._imageKey) {
            for (const mesh of this._imageMeshes) this.remove(mesh);
            this._imageMeshes = [];
            for (const [i, j, k] of rest) {
                const position = this._imagePosition([i, j, k]);
                for (const lobe of ['positive', 'negative']) {
                    const base = this.meshes?.[lobe];
                    if (!base) continue;
                    // Shared geometry AND material: the copies are the same
                    // surface seen in the next cell, so they must never drift
                    // from the original's colour, opacity or transparency
                    // policy — and sharing keeps a wide boundary cheap.
                    const image = new THREE.Mesh(base.geometry, base.material);
                    image.position.copy(position);
                    image.name = `${base.name}_image_${i}_${j}_${k}`;
                    image.userData.fieldLobe = lobe;
                    image.userData.isFieldPeriodicImage = true;
                    this._imageMeshes.push(image);
                    this.add(image);
                }
            }
            this._imageKey = key;
        }
        this._syncImageState();
    }

    /** Copies track their source mesh's geometry (replaced on every isovalue
     *  rebuild), visibility and render order. */
    _syncImageState() {
        for (const image of this._imageMeshes) {
            const base = this.meshes?.[image.userData.fieldLobe];
            if (!base) continue;
            image.geometry = base.geometry;
            image.material = base.material;
            image.visible = base.visible;
            image.renderOrder = base.renderOrder;
        }
    }

    /** Clip every copy to the boundary box (in the structure lattice's
     *  fractional coordinates). Nothing is clipped for the plain unit cell of
     *  the field's own grid — the field already ends at the cell faces there,
     *  so the default costs no clipping planes in the shader at all. With a
     *  different structure lattice even the unit cell clips, because the grid
     *  may reach past it (a cell reduced after the field was loaded). */
    _applyBoundsClipping(bounds) {
        const planes = isUnitBounds(bounds) && this._latticeIsOwnCell()
            ? null
            : makeFractionalBoundsClippingPlanes(this._boundaryVectors(), bounds);
        for (const key of ['positive', 'negative']) {
            const material = this.meshes?.[key]?.material;
            if (!material) continue;
            const before = material.clippingPlanes?.length ?? 0;
            material.clippingPlanes = planes;
            // The plane COUNT is compiled into the program.
            if (before !== (planes?.length ?? 0)) material.needsUpdate = true;
        }
    }

    /**
     * Point the marching-cubes backend at the field's current colour-by
     * field (or none), so the next pass interpolates it onto the vertices.
     */
    _syncColorField() {
        const colorBy = resolveColorBy(this.field);
        const active = this.marchingCubes.setColorField(colorBy?.field ?? null);
        this._colorSourceField = active ? colorBy.field : null;
    }

    /**
     * Colour both lobes from their stored per-vertex values and the field's
     * colorBy (colormap + range), or return them to the flat lobe colours
     * when colour-by is off. Cheap: no marching-cubes pass.
     *
     * @returns {boolean} false when colour-by is set but the geometry holds
     *   no values for that colour field (it changed since the last rebuild,
     *   or the backend cannot interpolate): the caller must rebuild.
     */
    applyColorBy() {
        const colorBy = resolveColorBy(this.field);
        let current = true;
        for (const key of ['positive', 'negative']) {
            const mesh = this.meshes?.[key];
            const geometry = mesh?.geometry;
            const material = mesh?.material;
            if (!geometry || !material) continue;
            const count = geometry.getAttribute('position')?.count ?? 0;
            const values = geometry.userData.colorByValues;
            let rgb = null;
            if (colorBy && count > 0) {
                if (values && values.length === count && geometry.userData.colorBySource === colorBy.field) {
                    rgb = mapValuesToColors(values, colorBy);
                } else {
                    current = false;
                }
            }
            geometry.userData.colorByRGB = rgb;
            // White under vertex colours so the colormap is not tinted; the
            // lobe colour otherwise (the live settings, see
            // setIsosurfaceMaterialSettings).
            material.color.copy(rgb ? _white : (key === 'positive' ? defaultPosColor : defaultNegColor));
            composeIsosurfaceVertexColors(mesh);
        }
        return current;
    }

    updateMesh(isoValue = this.field.isovalue, useAbsoluteIsoValue = false) {
        if (!groups.activeField) return;

        this._syncColorField();
        // Smoothing is global (see SMOOTHING_METHODS); optional so an older
        // wrapper without it still meshes plainly.
        this.marchingCubes?.setSmoothing?.(getIsosurfaceSmoothingSettings());
        let iso = isoValue;
        if (this.marchingCubes && this.meshes.positive && (iso >= 0 || useAbsoluteIsoValue)) {
            if (useAbsoluteIsoValue) {
                iso = Math.abs(isoValue);
            }
            this._lastIsoPositive = iso;
            this.marchingCubes.updateMesh(iso);
            const posData = this.marchingCubes.getVertices();
            //this.marchingCubes.sortVerticesToCamera(this.lastCameraPosition, posData.vertices, posData.normals);
            this.refreshGeometry("positive", posData.vertices, posData.normals, posData.values);
        }
        if (this.marchingCubes && this.meshes.negative && (iso < 0 || useAbsoluteIsoValue)) {
            if (useAbsoluteIsoValue) {
                iso = -Math.abs(isoValue);
            }
            this._lastIsoNegative = iso;
            this.marchingCubes.updateMesh(iso);
            const negData = this.marchingCubes.getVertices();
            //this.marchingCubes.sortVerticesToCamera(this.lastCameraPosition, negData.vertices, negData.normals);
            this.refreshGeometry("negative", negData.vertices, negData.normals, negData.values);
        }
        // Fresh geometry has no colour attribute yet: colour it (or confirm
        // the flat lobe colours) now. The focus alpha, if any, is layered on
        // by render/FocusRegionModule.js applyFocusToField afterwards.
        this.applyColorBy();
    }

    clearMesh() {
        // The copies share the geometry disposed below, so they go first; the
        // next setPeriodicBounds() rebuilds them against the new geometry.
        for (const image of this._imageMeshes) this.remove(image);
        this._imageMeshes = [];
        this._imageKey = null;
        if (this.meshes.positive) {
            this.remove(this.meshes.positive);
            this.meshes.positive.geometry.dispose();
        }
        if (this.meshes.negative) {
            this.remove(this.meshes.negative);
            this.meshes.negative.geometry.dispose();
        }
    }

    delete() {
        this.clearMesh();

        this.marchingCubes.delete();
    }

    setVisible(visible) {
        this.meshes.positive.visible = visible;
        this.meshes.negative.visible = visible;
        for (const image of this._imageMeshes) image.visible = visible;
        this.field.isVisible = visible;
    }
}