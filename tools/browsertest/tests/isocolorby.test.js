// Isosurface colour-by: an isosurface of field A coloured per vertex by the
// values of field B (Field.colorBy), the mechanism behind the NCI plot
// (s = 0.5 surface coloured by sign(λ2)ρ over ±0.04 with 'bgyor').
//
// Synthetic fields on YBCO's cell:
//   S — "s-like": 1 - gaussian, so the iso 0.5 surface is a sphere around the
//       LOW region (the NCI orientation: inside = below the isovalue);
//   L — signed, linear in x over [-0.1, 0.1] (the colour field);
//   M — on a different grid (must never be offered, and is ignored if set).
//
// Asserts:
//   (1) programmatic S.colorBy = {L, bgyor, ±0.04} + updateField gives an RGBA
//       vertex colour attribute with non-uniform RGB, blue at low x and red at
//       high x, on a white (untinted) material;
//   (2) the Field panel reflects a programmatic colour-by when S is selected
//       through the field list, and offers only same-grid fields;
//   (3) a colormap change through the dropdown remaps without a rebuild;
//   (4) None through the dropdown restores flat colours (attribute gone,
//       lobe colour back); picking L again through the dropdown re-colours
//       with the default symmetric range;
//   (5) with a focus region active the attribute carries colour-by RGB AND
//       focus alpha, whichever order they are applied in; clearing colour-by
//       keeps the focus alpha with white RGB;
//   (6) a mismatched-grid colour field is ignored.
'use strict';
const H = require('../harness');

(async () => {
  const { browser, page, errors } = await H.launchApp();
  H.check('webgl available', await H.webglAvailable(page));
  await H.loadDefaultStructure(page); // YBCO

  // --- setup: three fields, Field panel open, S selected ---------------------
  const setup = await page.evaluate(async () => {
    const { Field, FieldContainer } = await import('./model/index.js');
    const { fileBrowser, groups } = await import('./state/store.js');
    const { fieldBrowser } = await import('./ui/FieldPanel.js');
    const { setActiveField, updateField } = await import('./render/index.js');
    const { openPanel, refreshPanelAvailability } = await import('./ui/panels/PanelManager.js');

    const structure = fileBrowser.selectedStructure;
    const lat = structure.lattice;
    const make = (n, fill, label) => {
      const values = new Float32Array(n * n * n);
      let min = Infinity, max = -Infinity;
      for (let k = 0; k < n; k++)
        for (let j = 0; j < n; j++)
          for (let i = 0; i < n; i++) {
            const v = fill(i / (n - 1), j / (n - 1), k / (n - 1));
            values[i + n * (j + n * k)] = v;
            if (v < min) min = v;
            if (v > max) max = v;
          }
      const voxel = [0, 1, 2].map((axis) => lat[axis].map((c) => c / n));
      return new Field({
        nx: n, ny: n, nz: n, origin: [0, 0, 0], voxel, values, label,
        minValue: min, maxValue: max,
        absMinValue: Math.min(Math.abs(min), Math.abs(max)), absMaxValue: Math.max(Math.abs(min), Math.abs(max)),
        useAbsoluteIsoValue: false, isVisible: true,
      });
    };
    const S = make(24, (x, y, z) => 1 - Math.exp(-((x - 0.5) ** 2 + (y - 0.5) ** 2 + (z - 0.5) ** 2) / (2 * 0.18 ** 2)), 'S-like');
    const L = make(24, (x) => (x - 0.5) * 0.2, 'Linear colour');
    const M = make(12, (x) => x, 'Mismatched grid');
    S.isoValue = 0.5;
    L.isoValue = 0.01;
    M.isoValue = 0.5;
    structure.volumetricFields = new FieldContainer({ fileName: 'colorby.cube', source: 'Cube', fields: [S, L, M] });
    fieldBrowser.setCatalog(structure.volumetricFields.catalog);
    refreshPanelAvailability();
    openPanel('field');
    await new Promise((r) => setTimeout(r, 300));
    fieldBrowser.selectField(S);
    setActiveField(S, false);
    updateField(0.5);
    window.__colorBy = { S, L, M };

    const select = /** @type {HTMLSelectElement} */ (document.getElementById('fieldColorBySelect'));
    const mesh = groups.isosurfaceGroup.meshes.positive;
    return {
      hasSelect: !!select,
      options: select ? [...select.options].map((o) => o.textContent) : [],
      vertices: mesh.geometry.getAttribute('position')?.count ?? 0,
      hasColor: !!mesh.geometry.getAttribute('color'),
    };
  });
  H.check('setup: S builds a surface with flat colours and the colour-by dropdown exists',
    setup.hasSelect && setup.vertices > 0 && !setup.hasColor, JSON.stringify(setup));
  H.check('dropdown offers None + the same-grid field only (not itself, not the mismatched grid)',
    setup.options.length === 2 && /None/.test(setup.options[0]) && setup.options[1] === 'Linear colour',
    JSON.stringify(setup.options));

  // In-page probe of the live positive mesh's colours.
  const probe = () => page.evaluate(async () => {
    const { groups } = await import('./state/store.js');
    const { defaultPosColor } = await import('./model/Isosurface.js');
    const mesh = groups.isosurfaceGroup.meshes.positive;
    const geometry = mesh.geometry;
    const pos = geometry.getAttribute('position');
    const color = geometry.getAttribute('color');
    const out = {
      geometryId: geometry.id,
      vertices: pos.count,
      hasColor: !!color,
      itemSize: color?.itemSize ?? 0,
      vertexColors: mesh.material.vertexColors,
      materialHex: mesh.material.color.getHexString(),
      lobeHex: defaultPosColor.getHexString(),
      distinctRgb: 0, allRgbOne: null, minAlpha: null, maxAlpha: null,
      lowX: null, highX: null, meanRgb: null,
    };
    if (!color) return out;
    const a = color.array;
    const seen = new Set();
    let allOne = true, minA = Infinity, maxA = -Infinity;
    let lo = { x: Infinity, rgb: null }, hi = { x: -Infinity, rgb: null };
    const mean = [0, 0, 0];
    for (let i = 0; i < color.count; i++) {
      const r = a[i * 4], g = a[i * 4 + 1], b = a[i * 4 + 2], al = a[i * 4 + 3];
      seen.add(`${r.toFixed(2)},${g.toFixed(2)},${b.toFixed(2)}`);
      if (r !== 1 || g !== 1 || b !== 1) allOne = false;
      minA = Math.min(minA, al); maxA = Math.max(maxA, al);
      mean[0] += r; mean[1] += g; mean[2] += b;
      const x = pos.array[i * 3];
      if (x < lo.x) lo = { x, rgb: [r, g, b] };
      if (x > hi.x) hi = { x, rgb: [r, g, b] };
    }
    out.distinctRgb = seen.size;
    out.allRgbOne = allOne;
    out.minAlpha = minA;
    out.maxAlpha = maxA;
    out.lowX = lo.rgb;
    out.highX = hi.rgb;
    out.meanRgb = mean.map((v) => +(v / color.count).toFixed(3));
    return out;
  });

  // --- (1) programmatic colour-by -------------------------------------------
  await page.evaluate(async () => {
    const { updateField, requestRender } = await import('./render/index.js');
    const { S, L } = window.__colorBy;
    S.colorBy = { field: L, colormap: 'bgyor', min: -0.04, max: 0.04 };
    updateField();
    requestRender();
  });
  await page.waitForTimeout(400);
  const prog = await probe();
  H.check('programmatic colorBy: RGBA colour attribute with non-uniform RGB on a white material',
    prog.hasColor && prog.itemSize === 4 && prog.vertexColors === true && prog.distinctRgb > 10
      && prog.materialHex === 'ffffff' && prog.minAlpha === 1, JSON.stringify(prog));
  H.check('bgyor runs blue (negative, low x) to red (positive, high x)',
    prog.lowX && prog.highX && prog.lowX[2] > prog.lowX[0] && prog.highX[0] > prog.highX[2],
    JSON.stringify({ lowX: prog.lowX, highX: prog.highX }));
  await H.shotCanvas(page, 'isocolorby-bgyor');
  // Opaque copy for eyeballing the lighting: the s-like surface encloses the
  // LOW region, so its gradient normals point outward from the low side.
  await page.evaluate(async () => {
    const { setIsosurfaceMaterialSettings, applyMaterialSettingsToStoredIsosurfaces } = await import('./model/index.js');
    const { groups } = await import('./state/store.js');
    const { requestRender } = await import('./render/index.js');
    setIsosurfaceMaterialSettings({ opacity: 1 });
    applyMaterialSettingsToStoredIsosurfaces(groups.isosurfaceGroup, { opacity: 1 });
    requestRender();
  });
  await page.waitForTimeout(400);
  await H.shotCanvas(page, 'isocolorby-bgyor-opaque');
  await page.evaluate(async () => {
    const { setIsosurfaceMaterialSettings, applyMaterialSettingsToStoredIsosurfaces } = await import('./model/index.js');
    const { groups } = await import('./state/store.js');
    setIsosurfaceMaterialSettings({ opacity: 0.6 });
    applyMaterialSettingsToStoredIsosurfaces(groups.isosurfaceGroup, { opacity: 0.6 });
  });

  // --- (2) the panel reflects a programmatic colour-by on selection ----------
  const reflect = await page.evaluate(async () => {
    const radios = () => [...document.querySelectorAll('#fieldCatalogMount .fc-leaf-radio')];
    radios()[1].click(); // L
    await new Promise((r) => setTimeout(r, 200));
    const whileL = {
      value: /** @type {HTMLSelectElement} */ (document.getElementById('fieldColorBySelect')).value,
      optionsHidden: document.getElementById('fieldColorByOptions').hidden,
    };
    radios()[0].click(); // S, whose colorBy was set programmatically
    await new Promise((r) => setTimeout(r, 300));
    const select = /** @type {HTMLSelectElement} */ (document.getElementById('fieldColorBySelect'));
    return {
      whileL,
      value: select.value,
      selectedText: select.selectedOptions[0]?.textContent,
      optionsHidden: document.getElementById('fieldColorByOptions').hidden,
      colormap: /** @type {HTMLSelectElement} */ (document.getElementById('fieldColorByColormap')).value,
      bar: !!document.querySelector('#fieldColorByBarMount .cv-colorbar-wrapper'),
      barMin: /** @type {HTMLInputElement} */ (document.querySelector('#fieldColorByBarMount .cv-colorbar-value-input'))?.value,
    };
  });
  H.check('another field shows no colour-by (its own, unset)',
    reflect.whileL.value === '' && reflect.whileL.optionsHidden === true, JSON.stringify(reflect.whileL));
  H.check('selecting S through the list shows its programmatic colour-by, colormap and colour bar',
    reflect.selectedText === 'Linear colour' && reflect.optionsHidden === false && reflect.colormap === 'bgyor'
      && reflect.bar && Number(reflect.barMin) === -0.04, JSON.stringify(reflect));
  const afterReselect = await probe();
  H.check('re-selecting S rebuilds it still coloured', afterReselect.hasColor && afterReselect.distinctRgb > 10,
    JSON.stringify(afterReselect));

  // --- (3) colormap through the dropdown: remap, no rebuild -------------------
  await H.setSelect(page, 'fieldColorByColormap', 'viridis');
  await page.waitForTimeout(300);
  const viridis = await probe();
  H.check('colormap change remaps the same geometry (no marching-cubes rebuild)',
    viridis.geometryId === afterReselect.geometryId && viridis.hasColor
      && JSON.stringify(viridis.meanRgb) !== JSON.stringify(afterReselect.meanRgb),
    JSON.stringify({ before: afterReselect.meanRgb, after: viridis.meanRgb }));
  await H.setSelect(page, 'fieldColorByColormap', 'bgyor');

  // --- (4) None, then L again, through the dropdown ---------------------------
  await H.setSelect(page, 'fieldColorBySelect', '');
  await page.waitForTimeout(300);
  const none = await probe();
  const noneUi = await page.evaluate(() => ({
    optionsHidden: document.getElementById('fieldColorByOptions').hidden,
    bar: !!document.querySelector('#fieldColorByBarMount .cv-colorbar-wrapper'),
    colorBy: window.__colorBy.S.colorBy,
  }));
  H.check('None restores flat colours: no attribute, vertex colours off, lobe colour back',
    !none.hasColor && none.vertexColors === false && none.materialHex === none.lobeHex
      && noneUi.optionsHidden && !noneUi.bar && noneUi.colorBy === null,
    JSON.stringify({ none, noneUi }));
  await H.shotCanvas(page, 'isocolorby-flat');

  await H.setSelect(page, 'fieldColorBySelect', '0');
  await page.waitForTimeout(400);
  const again = await probe();
  const againState = await page.evaluate(() => {
    const cb = window.__colorBy.S.colorBy;
    return { field: cb?.field?.label, colormap: cb?.colormap, min: cb?.min, max: cb?.max };
  });
  H.check('picking L in the dropdown colours again with the symmetric default range',
    again.hasColor && again.distinctRgb > 10 && againState.field === 'Linear colour'
      && againState.min === -0.1 && againState.max === 0.1, JSON.stringify({ againState, again }));

  // --- (5) focus region coexistence -------------------------------------------
  const focusOn = await page.evaluate(async () => {
    const { fileBrowser } = await import('./state/store.js');
    const focus = await import('./render/FocusRegionModule.js');
    const structure = fileBrowser.selectedStructure;
    const wrapped = structure.periodic.visibleWrapped;
    const region = focus.createFocusRegion([{
      sourceIndex: wrapped.srcIndex[0], element: wrapped.elements[0], position: wrapped.cart[0],
    }]);
    region.innerRadius = 0.1;
    region.outerOpacity = 0.1;
    region.gradientEnabled = false;
    focus.applyFocusRegions();
  });
  void focusOn;
  await page.waitForTimeout(300);
  const focused = await probe();
  H.check('focus + colour-by: RGB from the colormap AND alpha from the focus',
    focused.hasColor && focused.distinctRgb > 10 && !focused.allRgbOne && focused.minAlpha <= 0.1001,
    JSON.stringify(focused));

  // Order: rebuild (colour-by first, then focus via updateField).
  await page.evaluate(async () => {
    const { updateField } = await import('./render/index.js');
    updateField();
  });
  const rebuilt = await probe();
  H.check('rebuild with focus active keeps both channels',
    rebuilt.hasColor && rebuilt.distinctRgb > 10 && rebuilt.minAlpha <= 0.1001, JSON.stringify(rebuilt));

  // Order: focus already present, colour-by removed then re-added (remap path).
  await H.setSelect(page, 'fieldColorBySelect', '');
  await page.waitForTimeout(300);
  const focusOnly = await probe();
  H.check('clearing colour-by under a focus keeps the alpha and turns RGB white',
    focusOnly.hasColor && focusOnly.allRgbOne === true && focusOnly.minAlpha <= 0.1001
      && focusOnly.materialHex === focusOnly.lobeHex, JSON.stringify(focusOnly));
  await H.setSelect(page, 'fieldColorBySelect', '0');
  await page.waitForTimeout(400);
  const both = await probe();
  H.check('re-adding colour-by under a focus restores the colormap and keeps the alpha',
    both.hasColor && both.distinctRgb > 10 && both.minAlpha <= 0.1001, JSON.stringify(both));
  await H.shotCanvas(page, 'isocolorby-focus');

  const focusOff = await page.evaluate(async () => {
    const focus = await import('./render/FocusRegionModule.js');
    focus.clearFocusRegions();
    focus.applyFocusRegions();
  });
  void focusOff;
  const noFocus = await probe();
  H.check('removing the focus keeps the colour-by attribute with alpha 1',
    noFocus.hasColor && noFocus.distinctRgb > 10 && noFocus.minAlpha === 1, JSON.stringify(noFocus));

  // --- tracer sanity: colour-by active under the ray tracer (flat colours there
  // in v1) must not break the encode or throw.
  await H.setSelect(page, 'renderPipelineMenu', 'raytrace');
  await page.waitForTimeout(1500);
  await H.shotCanvas(page, 'isocolorby-raytrace');
  const traced = await page.evaluate(async () => {
    const { general } = await import('./state/store.js');
    return { pipeline: general.renderPipeline };
  });
  H.check('ray tracer runs with a colour-by surface in the scene', traced.pipeline === 'raytrace',
    JSON.stringify(traced));
  await H.setSelect(page, 'renderPipelineMenu', 'depthpeel');
  await page.waitForTimeout(500);

  // --- (6) mismatched grid is ignored ------------------------------------------
  const mismatch = await page.evaluate(async () => {
    const { updateField } = await import('./render/index.js');
    const { groups } = await import('./state/store.js');
    const { S, M } = window.__colorBy;
    S.colorBy = { field: M, colormap: 'bgyor', min: 0, max: 1 };
    updateField();
    const mesh = groups.isosurfaceGroup.meshes.positive;
    const hasColor = !!mesh.geometry.getAttribute('color');
    S.colorBy = null;
    updateField();
    return { hasColor, vertexColors: mesh.material.vertexColors };
  });
  H.check('a colour field on a different grid is ignored (flat colours)',
    mismatch.hasColor === false && mismatch.vertexColors === false, JSON.stringify(mismatch));

  H.check('no page errors', errors.length === 0, errors.join(' | '));
  await H.finish(browser);
})().catch(H.crash);
