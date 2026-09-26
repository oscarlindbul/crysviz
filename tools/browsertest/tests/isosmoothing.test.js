// Isosurface smoothing (model/Isosurface.js SMOOTHING_METHODS /
// get/setIsosurfaceSmoothingSettings, applied by Isosurface.updateMesh(); the
// "Visual controls" section's #FieldSmoothingMethod select; ShareModule
// `smoothing`, v2.17).
//
// Synthetic fields on YBCO's cell, sampled so that grid point i sits at
// i/(n-1) — the normalised grid box the mesh positions live in — so analytic
// surfaces can be compared directly against the local vertex positions:
//   SPHERE — f = |p - c|, c = (0.5, 0.5, 0.5), iso R = 0.3: a closed sphere well
//            away from the cell faces;
//   SLAB   — f = z + wobble(x, y), iso 0.6: an open sheet crossing the x/y faces;
//   MASKED — the sphere with every point at x > 0.7 set above maskValue, so the
//            cubes touching them (x > 16/23) must stay empty.
//
// Asserts:
//   (1) the default method is 'off', the select shows "Off", SMOOTHING_METHODS
//       lists every method id, and 'off' reproduces the baseline mesh;
//   (2) every method on the sphere gives a sane non-indexed mesh (count > 0,
//       multiple of 3, position/normal same length, all finite); projected mesh
//       smoothers stay on the analytic sphere without widening its radial
//       spread; Taubin without projection does not shrink, Laplacian without
//       projection does; Loop / Catmull–Clark / tricubic multiply the triangle
//       count by ~4 / 6 / 4;
//   (3) Taubin keeps the rim vertices on the cell faces exactly where 'off'
//       puts them (no seams between periodic images) while moving the interior;
//   (4) gaussian and tricubic respect Field.maskValue (no vertices in masked
//       cubes, no crash);
//   (5) the UI select drives the settings, rebuilds the parameter controls and
//       re-meshes; Off restores the baseline vertex count;
//   (6) captureState/applySharedState round-trips the smoothing settings.
'use strict';
const H = require('../harness');

const N = 24;                    // grid points per axis
const R = 0.3;                   // sphere radius in the normalised grid box
const SPHERE_TOL = 0.015 * R;    // max |r - R| for projected mesh smoothers
const MASK_X = 16 / (N - 1);     // cubes beyond this x touch a masked point (i >= 17)

(async () => {
  const { browser, page, errors } = await H.launchApp();
  H.check('webgl available', await H.webglAvailable(page));
  await H.loadDefaultStructure(page); // YBCO

  // --- setup: in-page helpers, sphere field loaded, Field panel open ----------
  const setup = await page.evaluate(async ({ N, R, MASK_X }) => {
    const { Field, FieldContainer } = await import('./model/index.js');
    const iso = await import('./model/Isosurface.js');
    const { fileBrowser, groups } = await import('./state/store.js');
    const { fieldBrowser } = await import('./ui/FieldPanel.js');
    const { setActiveField, updateField } = await import('./render/index.js');
    const { openPanel, refreshPanelAvailability } = await import('./ui/panels/PanelManager.js');

    const structure = fileBrowser.selectedStructure;
    const lat = structure.lattice;
    const make = (fill, label, extra = {}) => {
      const values = new Float32Array(N * N * N);
      let min = Infinity, max = -Infinity;
      for (let k = 0; k < N; k++)
        for (let j = 0; j < N; j++)
          for (let i = 0; i < N; i++) {
            const v = fill(i / (N - 1), j / (N - 1), k / (N - 1));
            values[i + N * (j + N * k)] = v;
            if (v < min) min = v;
            if (v > max) max = v;
          }
      const voxel = [0, 1, 2].map((axis) => lat[axis].map((c) => c / N));
      return new Field({
        nx: N, ny: N, nz: N, origin: [0, 0, 0], voxel, values, label,
        minValue: min, maxValue: max,
        absMinValue: Math.min(Math.abs(min), Math.abs(max)), absMaxValue: Math.max(Math.abs(min), Math.abs(max)),
        useAbsoluteIsoValue: false, isVisible: true, ...extra,
      });
    };
    const dist = (x, y, z) => Math.hypot(x - 0.5, y - 0.5, z - 0.5);
    const sphere = make(dist, 'Sphere');
    sphere.isoValue = R;
    const slab = make((x, y, z) => 0.1 + z + 0.06 * Math.sin(2 * Math.PI * x) * Math.cos(2 * Math.PI * y), 'Slab');
    slab.isoValue = 0.6;
    const masked = make((x, y, z) => (x > 0.7 ? 10 : dist(x, y, z)), 'Masked sphere', { maskValue: 5 });
    masked.isoValue = R;
    if (masked.maskValue !== 5) masked.maskValue = 5;

    // Load one field as the active one through the same path as the field list.
    const load = async (field) => {
      structure.volumetricFields = new FieldContainer({ fileName: 'smoothing.cube', source: 'Cube', fields: [field] });
      fieldBrowser.setCatalog(structure.volumetricFields.catalog);
      refreshPanelAvailability();
      fieldBrowser.selectField(field);
      setActiveField(field, false);
      updateField(field.isoValue);
      await new Promise((r) => setTimeout(r, 100));
    };

    // Stats of the live positive mesh (local coords = normalised grid box).
    const stats = () => {
      const mesh = groups.isosurfaceGroup?.meshes?.positive;
      const g = mesh?.geometry;
      const pos = g?.getAttribute('position');
      const nor = g?.getAttribute('normal');
      const out = {
        count: pos?.count ?? 0,
        posLen: pos?.array?.length ?? 0,
        norLen: nor?.array?.length ?? -1,
        indexed: !!g?.index,
        finite: true,
        meanR: null, stdR: null, maxDev: null, maxX: null,
        boundary: [], interiorHash: 0,
      };
      if (!pos) return out;
      const p = pos.array;
      const n = nor?.array ?? [];
      for (let i = 0; i < p.length; i++) if (!Number.isFinite(p[i])) { out.finite = false; break; }
      for (let i = 0; i < n.length && out.finite; i++) if (!Number.isFinite(n[i])) { out.finite = false; break; }
      let s = 0, s2 = 0, maxDev = 0, maxX = -Infinity, hash = 0;
      const bset = new Set();
      const onFace = (v) => Math.abs(v) < 1e-6 || Math.abs(v - 1) < 1e-6;
      for (let v = 0; v < pos.count; v++) {
        const x = p[3 * v], y = p[3 * v + 1], z = p[3 * v + 2];
        const r = Math.hypot(x - 0.5, y - 0.5, z - 0.5);
        s += r; s2 += r * r;
        maxDev = Math.max(maxDev, Math.abs(r - R));
        maxX = Math.max(maxX, x);
        if (onFace(x) || onFace(y) || onFace(z)) bset.add(`${x.toFixed(5)},${y.toFixed(5)},${z.toFixed(5)}`);
        else hash += x * 1.3 + y * 2.7 + z * 3.1;
      }
      const m = s / pos.count;
      out.meanR = m;
      out.stdR = Math.sqrt(Math.max(0, s2 / pos.count - m * m));
      out.maxDev = maxDev;
      out.maxX = maxX;
      out.boundary = [...bset].sort();
      out.interiorHash = hash;
      return out;
    };

    // Apply a smoothing setting and rebuild the active field's surface.
    const apply = async (settings) => {
      iso.setIsosurfaceSmoothingSettings(settings);
      updateField(groups.activeField.isoValue);
      await new Promise((r) => setTimeout(r, 50));
      return stats();
    };

    window.__smooth = { iso, fields: { sphere, slab, masked }, load, stats, apply, MASK_X };

    // The field panel is only available once the structure has fields, so
    // load first, open the panel, then load again so it builds its controls.
    await load(sphere);
    openPanel('field');
    await new Promise((r) => setTimeout(r, 300));
    await load(sphere);

    const methods = iso.SMOOTHING_METHODS;
    const ids = Array.isArray(methods)
      ? methods.map((m) => (typeof m === 'string' ? m : m?.id))
      : Object.keys(methods ?? {});
    const select = /** @type {HTMLSelectElement} */ (document.getElementById('FieldSmoothingMethod'));
    return {
      ids,
      settings: iso.getIsosurfaceSmoothingSettings(),
      hasSelect: !!select,
      inMount: !!select?.closest('#fieldSmoothingMount'),
      inVisual: !!select?.closest('#fieldVisualContent'),
      visualToggle: !!document.getElementById('fieldVisualToggle'),
      selectValue: select?.value,
      selectText: select?.selectedOptions?.[0]?.textContent?.trim(),
      optionValues: select ? [...select.options].map((o) => o.value) : [],
      baseline: stats(),
    };
  }, { N, R, MASK_X });

  const ALL = ['off', 'laplacian', 'taubin', 'hc', 'loop', 'catmullClark', 'tricubic', 'gaussian'];
  const baseline = setup.baseline;
  H.check('SMOOTHING_METHODS lists every method id', ALL.every((id) => setup.ids.includes(id)),
    JSON.stringify(setup.ids));
  H.check('default smoothing method is off', setup.settings?.method === 'off', JSON.stringify(setup.settings));
  H.check('Visual controls section holds the smoothing select, showing "Off"',
    setup.hasSelect && setup.inMount && setup.inVisual && setup.visualToggle
      && setup.selectValue === 'off' && /^off$/i.test(setup.selectText || ''),
    JSON.stringify({ ...setup, baseline: undefined, settings: undefined, ids: undefined }));
  H.check('select offers every method', ALL.every((id) => setup.optionValues.includes(id)),
    JSON.stringify(setup.optionValues));
  H.check('baseline sphere: non-indexed mesh on the analytic sphere',
    baseline.count > 0 && baseline.count % 3 === 0 && !baseline.indexed && baseline.finite
      && baseline.posLen === baseline.norLen && baseline.maxDev < SPHERE_TOL && baseline.boundary.length === 0,
    JSON.stringify({ ...baseline, boundary: baseline.boundary.length }));

  // --- (1) explicit 'off' equals the baseline --------------------------------
  const off = await page.evaluate(() => window.__smooth.apply({ method: 'off' }));
  H.check("explicit 'off' reproduces the unsmoothed mesh",
    off.count === baseline.count && Math.abs(off.meanR - baseline.meanR) < 1e-9,
    JSON.stringify({ off: off.count, baseline: baseline.count }));

  // --- (2) every method on the sphere ----------------------------------------
  const CASES = {
    laplacian: { iterations: 10, lambda: 0.5, project: true },
    taubin: { iterations: 10, lambda: 0.5, passband: 0.1, project: true },
    hc: { iterations: 10, alpha: 0.1, beta: 0.5, project: true },
    loop: { levels: 1, project: true },
    catmullClark: { levels: 1, project: true },
    tricubic: { factor: 2, kernel: 'catmullRom' },
    gaussian: { sigma: 1 },
  };
  const res = {};
  for (const [method, params] of Object.entries(CASES)) {
    res[method] = await page.evaluate(({ method, params }) =>
      window.__smooth.apply({ method, params: { [method]: params } }), { method, params });
    const s = res[method];
    const settings = await page.evaluate(() => window.__smooth.iso.getIsosurfaceSmoothingSettings());
    H.check(`${method}: settings applied`, settings?.method === method, JSON.stringify(settings));
    H.check(`${method}: sane non-indexed mesh (count>0, x3, normals match, finite)`,
      s.count > 0 && s.count % 3 === 0 && !s.indexed && s.posLen === s.norLen && s.finite,
      JSON.stringify({ count: s.count, posLen: s.posLen, norLen: s.norLen, indexed: s.indexed, finite: s.finite }));
  }
  for (const method of ['laplacian', 'taubin', 'hc', 'loop', 'catmullClark']) {
    const s = res[method];
    H.check(`${method} (project): vertices stay on the analytic sphere`, s.maxDev < SPHERE_TOL,
      JSON.stringify({ maxDev: s.maxDev, tol: SPHERE_TOL }));
    // On an exact distance field plain marching cubes is already near-exact,
    // so only require the spread to stay small (projection lands on the
    // trilinear surface, slightly off the analytic one).
    H.check(`${method} (project): radial spread stays below 0.5% of R`,
      s.stdR <= 0.005 * R, JSON.stringify({ std: s.stdR, offStd: off.stdR }));
  }
  const triRatio = (s) => (s.count / 3) / (off.count / 3);
  H.check('loop levels 1: ~4x the triangles of off', Math.abs(triRatio(res.loop) - 4) <= 0.4,
    JSON.stringify({ ratio: triRatio(res.loop) }));
  H.check('catmullClark levels 1: ~6x the triangles of off', Math.abs(triRatio(res.catmullClark) - 6) <= 0.6,
    JSON.stringify({ ratio: triRatio(res.catmullClark) }));
  H.check('tricubic factor 2: roughly 4x the triangles of off',
    triRatio(res.tricubic) > 3 && triRatio(res.tricubic) < 5, JSON.stringify({ ratio: triRatio(res.tricubic) }));
  H.check('tricubic (catmullRom) keeps the sphere radius (mean within 2%)',
    Math.abs(res.tricubic.meanR - R) < 0.02 * R, JSON.stringify({ meanR: res.tricubic.meanR }));
  H.check('gaussian changes the surface but keeps it near the sphere (mean radius within 10%)',
    Math.abs(res.gaussian.meanR - R) < 0.1 * R, JSON.stringify({ meanR: res.gaussian.meanR }));

  const noProj = await page.evaluate(async () => {
    const { apply } = window.__smooth;
    const taubin = await apply({ method: 'taubin', params: { taubin: { iterations: 10, lambda: 0.5, passband: 0.1, project: false } } });
    const laplacian = await apply({ method: 'laplacian', params: { laplacian: { iterations: 50, lambda: 0.5, project: false } } });
    return { taubin: taubin.meanR, taubinCount: taubin.count, laplacian: laplacian.meanR };
  });
  H.check('taubin without projection does not shrink (mean radius within 2% of off)',
    noProj.taubinCount > 0 && Math.abs(noProj.taubin - off.meanR) < 0.02 * off.meanR,
    JSON.stringify({ ...noProj, off: off.meanR }));
  H.check('laplacian without projection (50 iterations) shrinks the sphere',
    noProj.laplacian < off.meanR * 0.99, JSON.stringify({ ...noProj, off: off.meanR }));
  await H.shotCanvas(page, 'isosmoothing-sphere');

  // --- (3) rim vertices on the cell faces stay put ----------------------------
  const rim = await page.evaluate(async () => {
    const { load, apply, fields } = window.__smooth;
    window.__smooth.iso.setIsosurfaceSmoothingSettings({ method: 'off' });
    await load(fields.slab);
    const plain = await apply({ method: 'off' });
    const smooth = await apply({ method: 'taubin', params: { taubin: { iterations: 10, lambda: 0.5, passband: 0.1, project: false } } });
    const b1 = new Set(plain.boundary);
    return {
      offCount: plain.count, count: smooth.count, finite: smooth.finite,
      offRim: plain.boundary.length, rim: smooth.boundary.length,
      same: plain.boundary.length === smooth.boundary.length && smooth.boundary.every((k) => b1.has(k)),
      missing: plain.boundary.filter((k) => !smooth.boundary.includes(k)).slice(0, 5),
      interiorMoved: Math.abs(plain.interiorHash - smooth.interiorHash) > 1e-4,
    };
  });
  H.check('slab: the surface reaches the cell faces', rim.offCount > 0 && rim.offRim > 0, JSON.stringify(rim));
  H.check('taubin keeps the cell-face rim vertices exactly where off puts them',
    rim.count > 0 && rim.finite && rim.same, JSON.stringify(rim));
  H.check('taubin still moves the interior of the open surface', rim.interiorMoved, JSON.stringify(rim));

  // --- (4) masked field -------------------------------------------------------
  const mask = await page.evaluate(async () => {
    const { load, apply, fields, MASK_X } = window.__smooth;
    window.__smooth.iso.setIsosurfaceSmoothingSettings({ method: 'off' });
    await load(fields.masked);
    const out = {};
    const pick = (s) => ({ count: s.count, finite: s.finite, maxX: s.maxX, x3: s.count % 3 === 0 });
    out.off = pick(await apply({ method: 'off' }));
    out.gaussian = pick(await apply({ method: 'gaussian', params: { gaussian: { sigma: 1 } } }));
    out.tricubic = pick(await apply({ method: 'tricubic', params: { tricubic: { factor: 2, kernel: 'catmullRom' } } }));
    out.tricubicB = pick(await apply({ method: 'tricubic', params: { tricubic: { factor: 2, kernel: 'bspline' } } }));
    out.limit = MASK_X;
    return out;
  });
  H.check('masked field (off): no vertices in masked cubes', mask.off.count > 0 && mask.off.maxX <= mask.limit + 1e-4,
    JSON.stringify(mask.off));
  for (const key of ['gaussian', 'tricubic', 'tricubicB']) {
    const s = mask[key];
    H.check(`masked field (${key}): mesh built, no vertices in masked cubes`,
      s.count > 0 && s.x3 && s.finite && s.maxX <= mask.limit + 1e-4, JSON.stringify({ ...s, limit: mask.limit }));
  }

  // --- (5) UI -----------------------------------------------------------------
  await page.evaluate(async () => {
    const { load, fields, iso } = window.__smooth;
    iso.setIsosurfaceSmoothingSettings({ method: 'off' });
    await load(fields.sphere);
    // Expand Visual controls if it is collapsed (the select works either way).
    const content = document.getElementById('fieldVisualContent');
    const hidden = content && (content.hidden || getComputedStyle(content).display === 'none');
    if (hidden) document.getElementById('fieldVisualToggle')?.click();
  });
  await page.waitForTimeout(200);
  const beforeUi = await page.evaluate(() => {
    const mount = document.getElementById('fieldSmoothingMount');
    return {
      offControls: mount ? mount.querySelectorAll('input, select:not(#FieldSmoothingMethod)').length : -1,
      count: window.__smooth.stats().count,
    };
  });
  await H.setSelect(page, 'FieldSmoothingMethod', 'loop');
  await page.waitForTimeout(500);
  const uiLoop = await page.evaluate(() => {
    const mount = document.getElementById('fieldSmoothingMount');
    return {
      method: window.__smooth.iso.getIsosurfaceSmoothingSettings().method,
      controls: mount ? mount.querySelectorAll('input, select:not(#FieldSmoothingMethod)').length : 0,
      count: window.__smooth.stats().count,
    };
  });
  H.check('UI: picking Loop sets the method and shows parameter controls',
    uiLoop.method === 'loop' && uiLoop.controls > 0, JSON.stringify({ beforeUi, uiLoop }));
  H.check('UI: picking Loop re-meshes (~4x vertices)',
    Math.abs(uiLoop.count / beforeUi.count - 4) <= 0.4, JSON.stringify({ beforeUi, uiLoop }));
  await H.setSelect(page, 'FieldSmoothingMethod', 'tricubic');
  await page.waitForTimeout(500);
  const uiTri = await page.evaluate(() => {
    const mount = document.getElementById('fieldSmoothingMount');
    return {
      method: window.__smooth.iso.getIsosurfaceSmoothingSettings().method,
      controls: mount ? mount.querySelectorAll('input, select:not(#FieldSmoothingMethod)').length : 0,
    };
  });
  H.check('UI: picking Tricubic swaps in its own parameter controls',
    uiTri.method === 'tricubic' && uiTri.controls > 0, JSON.stringify(uiTri));
  await H.setSelect(page, 'FieldSmoothingMethod', 'off');
  await page.waitForTimeout(500);
  const uiOff = await page.evaluate(() => ({
    method: window.__smooth.iso.getIsosurfaceSmoothingSettings().method,
    count: window.__smooth.stats().count,
  }));
  H.check('UI: switching back to Off restores the baseline vertex count',
    uiOff.method === 'off' && uiOff.count === baseline.count, JSON.stringify({ uiOff, baseline: baseline.count }));

  // --- (6) save / restore -----------------------------------------------------
  const rt = await page.evaluate(async () => {
    const { iso, stats } = window.__smooth;
    const { captureState, applySharedState } = await import('./ui/ShareModule.js');
    const { updateField } = await import('./render/index.js');
    const { groups } = await import('./state/store.js');
    iso.setIsosurfaceSmoothingSettings({
      method: 'taubin', params: { taubin: { iterations: 7, lambda: 0.42, passband: 0.13, project: false } },
    });
    updateField(groups.activeField.isoValue);
    const state = captureState({ includeFrames: true, includeFields: true });
    const copy = JSON.parse(JSON.stringify(state));
    const saved = copy.smoothing ?? copy.fields?.smoothing ?? null;

    iso.setIsosurfaceSmoothingSettings({ method: 'off' });
    applySharedState(copy, 'smoothing.crysviz');
    await new Promise((r) => setTimeout(r, 500));
    // Loading a state rebuilds the panels; reopen the field panel, whose
    // control is built from the (restored) global settings.
    const { openPanel } = await import('./ui/panels/PanelManager.js');
    openPanel('field');
    await new Promise((r) => setTimeout(r, 300));
    const restored = iso.getIsosurfaceSmoothingSettings();
    return {
      version: state.version,
      saved,
      restored,
      selectValue: /** @type {HTMLSelectElement} */ (document.getElementById('FieldSmoothingMethod'))?.value,
      count: stats().count,
    };
  });
  const [maj, min] = String(rt.version || '0.0').split('.').map(Number);
  H.check('captureState version is at least 2.17', maj > 2 || (maj === 2 && min >= 17), rt.version);
  H.check('captureState carries the smoothing settings',
    rt.saved?.method === 'taubin' && rt.saved?.params?.taubin?.iterations === 7, JSON.stringify(rt.saved));
  const tp = rt.restored?.params?.taubin;
  H.check('applySharedState restores method + params',
    rt.restored?.method === 'taubin' && tp?.iterations === 7 && Math.abs(tp?.lambda - 0.42) < 1e-9
      && Math.abs(tp?.passband - 0.13) < 1e-9 && tp?.project === false, JSON.stringify(rt.restored));
  H.check('restore updates the select and rebuilds a surface', rt.selectValue === 'taubin' && rt.count > 0,
    JSON.stringify({ selectValue: rt.selectValue, count: rt.count }));

  // --- teardown: back to off --------------------------------------------------
  const reset = await page.evaluate(async () => {
    const { updateField } = await import('./render/index.js');
    const { groups } = await import('./state/store.js');
    const { iso } = window.__smooth;
    iso.setIsosurfaceSmoothingSettings({ method: 'off' });
    if (groups.activeField && groups.isosurfaceGroup) updateField(groups.activeField.isoValue);
    return iso.getIsosurfaceSmoothingSettings().method;
  });
  H.check('reset to off', reset === 'off', reset);

  H.check('no page errors', errors.length === 0, errors.join(' | '));
  await H.finish(browser);
})().catch(H.crash);
