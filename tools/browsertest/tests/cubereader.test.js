// Gaussian cube files through the real loader (core/crystal-viewer.js
// loadStructure → formats.js → io/ReadCubeModule.js → the field panel and the
// marching-cubes isosurface).
//
// The main case is a synthetic Gaussian blob written with NEGATIVE grid counts
// (Å units) and a `.cub` name — the two things the old reader got wrong (it
// multiplied the axes by the signed count, flipping the cell) or did not
// recognise. The node unit tests (tools/unittest/cubereader.test.mjs) pin the
// parser in detail; this checks that its output drives the app: the field
// label from line 2, the atoms, a right-handed cell, the new Field unit and
// periodicity flags, and an isosurface that actually renders. Two of the unit
// fixtures (a Gradient cube and an orbital cube) check the multi-field paths.
//
// Periodicity follows the feature branch's meaning: a cube loads periodic by
// default; `buildCubeStructure(..., { periodic: false })` (the block prompt,
// or `loadStructure(..., { periodic: false })`) makes every field a finite
// block, MO and NVal = 4 cubes included. A .crysviz save and a share payload
// carry `valueUnit` and `periodic` back.
'use strict';
const fs = require('fs');
const path = require('path');
const H = require('../harness');
const { PNG } = require(`${__dirname}/../env/node_modules/pngjs`);

/** Pixels that differ substantially between two screenshots. */
function changedPixelCount(fileA, fileB) {
  const a = PNG.sync.read(fs.readFileSync(fileA));
  const b = PNG.sync.read(fs.readFileSync(fileB));
  let n = 0;
  const total = Math.min(a.width * a.height, b.width * b.height);
  for (let i = 0; i < total; i++) {
    const o = i * 4;
    const d = Math.abs(a.data[o] - b.data[o]) + Math.abs(a.data[o + 1] - b.data[o + 1])
      + Math.abs(a.data[o + 2] - b.data[o + 2]);
    if (d > 90) n++;
  }
  return n;
}

const UNIT_FIXTURES = path.join(__dirname, '..', '..', 'unittest', 'fixtures', 'cube');
const fixture = (name) => fs.readFileSync(path.join(UNIT_FIXTURES, name), 'utf8');

/** A 24³ Gaussian blob in a 6 Å box, Å units (negative counts), 2 atoms. */
function blobCube() {
  const n = 24;
  const step = 0.25;
  const f = (x) => x.toFixed(6).padStart(12);
  const e = (x) => x.toExponential(5).toUpperCase().padStart(13);
  const lines = [
    ' Synthetic blob',
    ' Electron density from Total SCF Density',
    `    2${f(-3)}${f(-3)}${f(-3)}    1`,
    `  -${n}${f(step)}${f(0)}${f(0)}`,
    `  -${n}${f(0)}${f(step)}${f(0)}`,
    `  -${n}${f(0)}${f(0)}${f(step)}`,
    `    8${f(8)}${f(0)}${f(0)}${f(0)}`,
    `    1${f(1)}${f(0.9)}${f(0)}${f(0)}`,
  ];
  const c = (n - 1) / 2;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const rec = [];
      for (let k = 0; k < n; k++) {
        const r2 = ((i - c) ** 2 + (j - c) ** 2 + (k - c) ** 2) * step * step;
        rec.push(Math.exp(-r2 / 2));
      }
      for (let s = 0; s < rec.length; s += 6) lines.push(rec.slice(s, s + 6).map(e).join(''));
    }
  }
  return lines.join('\n') + '\n';
}

(async () => {
  const { browser, page, errors } = await H.launchApp();
  H.check('webgl available', await H.webglAvailable(page));

  const load = (text, name) => page.evaluate(async ({ text, name }) => {
    const cv = await import('./core/crystal-viewer.js');
    const { fileBrowser, groups } = await import('./state/store.js');
    const { detectFormat, headOf } = await import('./io/index.js');
    await cv.loadStructure(text, name);
    const s = fileBrowser.selectedStructure;
    const vf = s?.volumetricFields;
    const iso = groups.isosurfaceGroup;
    return {
      detected: detectFormat({ fileName: name, head: headOf(text) }).id,
      atoms: s?.atoms?.length ?? 0,
      elements: s?.elements ?? [],
      lattice: s?.lattice ?? null,
      source: vf?.source ?? null,
      fields: (vf?.fields ?? []).map((f) => ({
        label: f.label, valueUnit: f.valueUnit, periodic: f.periodic, n: [f.nx, f.ny, f.nz],
        max: f.maxValue, voxel: f.voxel,
      })),
      activeLabel: groups.activeField?.label ?? null,
      posVerts: iso?.meshes?.positive?.geometry?.attributes?.position?.count ?? 0,
      inScene: !!iso?.parent,
    };
  }, { text, name });

  // --- the blob: Å units, .cub name ------------------------------------------
  const blob = await load(blobCube(), 'blob.cub');
  await page.waitForTimeout(1000);
  H.check('a .cub file is detected as a cube', blob.detected === 'cube', blob.detected);
  H.check('two atoms, O and H', blob.atoms === 2 && blob.elements.join(',') === 'O,H', JSON.stringify(blob.elements));
  const lat = blob.lattice || [];
  H.check('Å axes are not flipped: lattice = +6 Å along each axis',
    Math.abs(lat[0]?.[0] - 6) < 1e-6 && Math.abs(lat[1]?.[1] - 6) < 1e-6 && Math.abs(lat[2]?.[2] - 6) < 1e-6,
    JSON.stringify(lat));
  const [f0] = blob.fields;
  H.check('one field, labelled from line 2', blob.fields.length === 1
    && f0.label === 'Electron density from Total SCF Density', JSON.stringify(blob.fields.map((f) => f.label)));
  H.check('cube field carries valueUnit e/bohr³ and loads periodic by default',
    f0?.valueUnit === 'e/bohr³' && f0?.periodic === true, JSON.stringify(f0));
  H.check('container source is Cube', blob.source === 'Cube', blob.source);
  H.check('the field is active and its isosurface has vertices',
    blob.activeLabel === f0?.label && blob.posVerts > 0 && blob.inScene, JSON.stringify({ a: blob.activeLabel, v: blob.posVerts }));
  // Pixels, not just a mesh: the same view with the isosurface group hidden
  // must differ substantially from the one with it shown.
  const setIsoVisible = (visible) => page.evaluate(async (visible) => {
    const { groups } = await import('./state/store.js');
    const { requestRender } = await import('./render/index.js');
    if (groups.isosurfaceGroup) groups.isosurfaceGroup.visible = visible;
    requestRender();
  }, visible);
  const withIso = await H.shotCanvas(page, 'cube-blob');
  await setIsoVisible(false);
  await page.waitForTimeout(800);
  const withoutIso = await H.shotCanvas(page, 'cube-blob-hidden');
  await setIsoVisible(true);
  const delta = changedPixelCount(withIso, withoutIso);
  H.check('the blob isosurface is drawn', delta > 500, JSON.stringify({ delta }));

  // --- Gradient cube: four fields -------------------------------------------
  const grad = await load(fixture('gradient_nval4.cube'), 'grad.cube');
  H.check('NVal = 4 gradient cube gives Density, ∂ρ/∂x, ∂ρ/∂y, ∂ρ/∂z',
    JSON.stringify(grad.fields.map((f) => f.label)) === JSON.stringify(['Density', '∂ρ/∂x', '∂ρ/∂y', '∂ρ/∂z']),
    JSON.stringify(grad.fields.map((f) => f.label)));
  H.check('only the density of a gradient cube has a unit',
    JSON.stringify(grad.fields.map((f) => f.valueUnit)) === JSON.stringify(['e/bohr³', null, null, null]),
    JSON.stringify(grad.fields.map((f) => f.valueUnit)));

  // --- orbital cube: one field per MO ---------------------------------------
  const mo = await load(fixture('mo_two_orbitals.cube'), 'orbitals.cube');
  H.check('orbital cube gives MO 5 and MO 6', JSON.stringify(mo.fields.map((f) => f.label)) === '["MO 5","MO 6"]',
    JSON.stringify(mo.fields.map((f) => f.label)));
  H.check('orbital cube container says so', mo.source === 'Cube (orbitals)', mo.source);

  // --- periodic vs block on every field (buildCubeStructure / readCubeFile) --
  const flags = await page.evaluate(async (texts) => {
    const { buildCubeStructure, readCubeFile } = await import('./io/ReadCubeModule.js');
    const { parseCube } = await import('./io/cubeParse.js');
    const out = {};
    for (const [name, text] of Object.entries(texts)) {
      const summary = (r) => {
        const fs = r.structure_with_field.volumetricFields.fields;
        return {
          periodic: fs.map((f) => f.periodic),
          origins: fs.map((f) => JSON.stringify(f.origin)),
          labels: fs.map((f) => f.label),
          units: fs.map((f) => f.valueUnit),
          source: r.structure_with_field.volumetricFields.source,
        };
      };
      out[name] = {
        read: summary(readCubeFile(text, name)),
        periodic: summary(buildCubeStructure(parseCube(text), name, { periodic: true })),
        block: summary(buildCubeStructure(parseCube(text), name, { periodic: false })),
      };
    }
    return out;
  }, { blob: blobCube(), grad: fixture('gradient_nval4.cube'), mo: fixture('mo_two_orbitals.cube') });
  for (const [name, r] of Object.entries(flags)) {
    H.check(`${name}: readCubeFile and { periodic: true } give periodic === true on every field`,
      [...r.read.periodic, ...r.periodic.periodic].every((p) => p === true), JSON.stringify(r));
    H.check(`${name}: { periodic: false } gives periodic === false on every field, one shared origin`,
      r.block.periodic.length > 0 && r.block.periodic.every((p) => p === false)
        && new Set(r.block.origins).size === 1, JSON.stringify(r.block));
    H.check(`${name}: labels and units do not depend on the layout`,
      JSON.stringify(r.block.labels) === JSON.stringify(r.read.labels)
        && JSON.stringify(r.block.units) === JSON.stringify(r.read.units), JSON.stringify(r));
  }
  H.check('MO cube built as a block keeps one field per orbital and its source',
    flags.mo.block.periodic.length === 2 && flags.mo.block.source === 'Cube (orbitals)', JSON.stringify(flags.mo.block));
  H.check('NVal = 4 cube built as a block keeps four fields', flags.grad.block.periodic.length === 4,
    JSON.stringify(flags.grad.block));

  // --- .crysviz save -> load, and the share payload restore -----------------
  const roundTrip = await page.evaluate(async (text) => {
    const cv = await import('./core/crystal-viewer.js');
    const { captureState, applySharedState, waitForStateRestoration } = await import('./ui/ShareModule.js');
    const { fileBrowser } = await import('./state/store.js');
    const fieldsNow = () => fileBrowser.selectedStructure.volumetricFields.fields
      .map((f) => ({ periodic: f.periodic, valueUnit: f.valueUnit, origin: f.origin }));
    const out = {};
    for (const periodic of [true, false]) {
      await cv.loadStructure(text, `rt-${periodic}.cub`, false, '', { periodic });
      const before = fieldsNow();
      // What the Save panel's .crysviz button writes (ui/SavePanel.js).
      const saved = JSON.stringify({ format: 'crysviz', ...captureState({ includeFrames: true, includeFields: true }) }, null, 2);
      const savedFields = JSON.parse(saved).fields.fields;
      // Each key once per field in the saved text, not only in the parsed object.
      const keyCounts = ['"periodic": ', '"valueUnit": '].map((k) => saved.split(k).length - 1);
      await cv.loadStructure(saved, `rt-${periodic}.crysviz`);
      out[periodic] = { before, savedFields: savedFields.map((f) => ({ periodic: f.periodic, valueUnit: f.valueUnit })),
        keyCounts, fieldCount: savedFields.length, after: fieldsNow() };
    }
    // The share-payload restorer with tampered values (share links are untrusted).
    const state = captureState({ includeFrames: true, includeFields: true });
    const tampered = {};
    for (const [key, value] of Object.entries({ str: 'false', zero: 0, literal: false })) {
      const copy = JSON.parse(JSON.stringify(state));
      copy.fields.fields[0].periodic = value;
      if (key === 'str') copy.fields.fields[0].valueUnit = { evil: true };
      applySharedState(copy);
      await waitForStateRestoration();
      tampered[key] = fieldsNow()[0];
    }
    return { out, tampered };
  }, blobCube());
  for (const periodic of ['true', 'false']) {
    const r = roundTrip.out[periodic];
    H.check(`.crysviz round trip (periodic ${periodic}) keeps valueUnit and periodic`,
      r.after.length === r.before.length && r.after.every((f, i) => f.periodic === r.before[i].periodic
        && f.valueUnit === r.before[i].valueUnit && f.valueUnit === 'e/bohr³'
        && f.periodic === (periodic === 'true')), JSON.stringify(r));
    H.check(`.crysviz (periodic ${periodic}) emits periodic and valueUnit once per field`,
      r.keyCounts.every((n) => n === r.fieldCount), JSON.stringify(r.keyCounts));
  }
  H.check('share payload: periodic "false" (a string) or 0 restores a periodic field',
    roundTrip.tampered.str.periodic === true && roundTrip.tampered.zero.periodic === true,
    JSON.stringify(roundTrip.tampered));
  H.check('share payload: a literal false restores a block', roundTrip.tampered.literal.periodic === false,
    JSON.stringify(roundTrip.tampered.literal));
  H.check('share payload: a non-string valueUnit is dropped', roundTrip.tampered.str.valueUnit === null,
    JSON.stringify(roundTrip.tampered.str));
  H.check('share payload: a string valueUnit survives', roundTrip.tampered.zero.valueUnit === 'e/bohr³',
    JSON.stringify(roundTrip.tampered.zero));

  // --- a truncated file is refused, not half-loaded -------------------------
  const truncated = await page.evaluate(async (text) => {
    const cv = await import('./core/crystal-viewer.js');
    try {
      await cv.loadStructure(text, 'truncated.cube');
      return 'loaded';
    } catch (e) {
      return String(e?.message || e);
    }
  }, fixture('truncated.cube'));
  H.check('a truncated cube throws a clear error', /truncated/.test(truncated), truncated);
  // loadStructure reports load failures on the console too; those are expected here.
  const unexpected = errors.filter((e) => !/truncated/.test(e));

  H.check('no page errors', unexpected.length === 0, unexpected.join('\n'));
  await H.finish(browser);
})().catch(H.crash);
