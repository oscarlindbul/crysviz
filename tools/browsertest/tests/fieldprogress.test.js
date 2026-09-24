// The progress bar under the volumetric field list (ui/FieldProgressWidget.js
// over the task registry in state/fieldTasks.js).
//
// Checks: the bar sits directly under the list; a task is shown once it
// outlasts the show delay (at once when `immediate`), carries its label and
// stage, moves, and fades away when done; a quick task never shows; the real
// producers report tasks (SCF-NCI, "Combine fields" on a large grid, loading a
// cube file while the panel is open); and a panel rebuilt mid-task shows the
// task again.
'use strict';
const H = require('../harness');

/** A small Gaussian blob cube (bohr), one O atom, so NCI has a density. */
function blobCube(n = 24, name = 'Blob') {
  const step = 0.3;
  const f = (x) => x.toFixed(6).padStart(12);
  const e = (x) => x.toExponential(5).toUpperCase().padStart(13);
  const lines = [
    ` ${name}`,
    ' Electron density from Total SCF Density',
    `    1${f(0)}${f(0)}${f(0)}    1`,
    `  ${n}${f(step)}${f(0)}${f(0)}`,
    `  ${n}${f(0)}${f(step)}${f(0)}`,
    `  ${n}${f(0)}${f(0)}${f(step)}`,
    `    8${f(8)}${f((n - 1) * step / 2)}${f((n - 1) * step / 2)}${f((n - 1) * step / 2)}`,
  ];
  const c = (n - 1) / 2;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const rec = [];
      for (let k = 0; k < n; k++) {
        const r2 = ((i - c) ** 2 + (j - c) ** 2 + (k - c) ** 2) * step * step;
        rec.push(0.3 * Math.exp(-r2 / 3));
      }
      for (let s = 0; s < rec.length; s += 6) lines.push(rec.slice(s, s + 6).map(e).join(''));
    }
  }
  return lines.join('\n') + '\n';
}

/** Current state of the progress block. */
const readBar = (page) => page.evaluate(() => {
  const root = document.querySelector('#fieldCatalogMount .field-progress');
  const rows = [...(root?.querySelectorAll('.field-progress-row') || [])];
  return {
    present: Boolean(root),
    hidden: root ? root.hidden : true,
    underList: Boolean(root && root.previousElementSibling?.classList.contains('field-catalog-tree')),
    rows: rows.map((r) => ({
      text: r.querySelector('.field-progress-label')?.textContent || '',
      value: Number(r.querySelector('progress')?.value),
      done: r.classList.contains('field-progress-row-done'),
    })),
  };
});

(async () => {
  const { browser, page, errors } = await H.launchApp();
  H.check('webgl available', await H.webglAvailable(page));

  await page.evaluate(async (text) => {
    const cv = await import('./core/crystal-viewer.js');
    await cv.loadStructure(text, 'blob.cube');
    const { openPanel } = await import('./ui/panels/PanelManager.js');
    openPanel('field');
    // Every task any producer reports, for the checks further down.
    const { subscribeFieldTasks, fieldTasks } = await import('./state/fieldTasks.js');
    window.__taskLog = [];
    subscribeFieldTasks(() => {
      for (const t of fieldTasks()) window.__taskLog.push(`${t.label}|${t.state}|${t.immediate}`);
    });
    // Rows as they are added, so a short-lived one is still caught.
    window.__rowsSeen = [];
    new MutationObserver((records) => {
      for (const r of records) {
        for (const node of r.addedNodes) {
          if (node.classList?.contains('field-progress-row')) {
            // The label is filled right after the row is appended.
            queueMicrotask(() => window.__rowsSeen.push(
              node.querySelector('.field-progress-label')?.textContent || ''));
          }
        }
      }
    }).observe(document.body, { childList: true, subtree: true });
  }, blobCube());
  await page.waitForTimeout(2000);

  const idle = await readBar(page);
  H.check('the progress block sits directly under the field list', idle.present && idle.underList,
    JSON.stringify(idle));
  H.check('with nothing running it is hidden', idle.hidden && idle.rows.length === 0, JSON.stringify(idle));

  // --- an immediate task: shown at once, moves, fades, goes -----------------
  await page.evaluate(async () => {
    const { startFieldTask } = await import('./state/fieldTasks.js');
    window.__t = startFieldTask('Test task', { immediate: true });
    window.__t.setStage('stage one', 0, 0.5);
  });
  await page.waitForTimeout(150);
  const a1 = await readBar(page);
  await page.waitForTimeout(600);
  const a2 = await readBar(page);
  H.check('an immediate task shows its label and stage', !a1.hidden && a1.rows.length === 1
    && a1.rows[0].text === 'Test task — stage one', JSON.stringify(a1));
  H.check('the bar creeps forward inside the stage range',
    a2.rows[0]?.value > a1.rows[0]?.value && a2.rows[0]?.value < 0.5, JSON.stringify([a1.rows, a2.rows]));
  await page.evaluate(() => window.__t.setProgress(0.7));
  await page.waitForTimeout(150);
  const a3 = await readBar(page);
  H.check('reported progress sets the bar', Math.abs(a3.rows[0]?.value - 0.7) < 1e-6, JSON.stringify(a3.rows));
  await page.evaluate(() => window.__t.done());
  await page.waitForTimeout(50);
  const a4 = await readBar(page);
  H.check('a finished task fills the bar and fades', a4.rows[0]?.value === 1 && a4.rows[0]?.done,
    JSON.stringify(a4.rows));
  await page.waitForTimeout(800);
  const a5 = await readBar(page);
  H.check('then it disappears and the block hides again', a5.hidden && a5.rows.length === 0, JSON.stringify(a5));

  // --- the show delay ---------------------------------------------------------
  const quick = await page.evaluate(async () => {
    const { startFieldTask } = await import('./state/fieldTasks.js');
    const t = startFieldTask('Quick task');
    const seen = [];
    for (const wait of [50, 100]) {
      await new Promise((r) => { setTimeout(r, wait); });
      seen.push(!document.querySelector('#fieldCatalogMount .field-progress').hidden);
    }
    t.done();
    await new Promise((r) => { setTimeout(r, 400); });
    seen.push(!document.querySelector('#fieldCatalogMount .field-progress').hidden);
    return seen;
  });
  H.check('a task that finishes inside the show delay never shows', quick.every((v) => !v), JSON.stringify(quick));

  const slow = await page.evaluate(async () => {
    const { startFieldTask } = await import('./state/fieldTasks.js');
    const t = startFieldTask('Slow task');
    const visible = () => !document.querySelector('#fieldCatalogMount .field-progress').hidden;
    await new Promise((r) => { setTimeout(r, 100); });
    const early = visible();
    await new Promise((r) => { setTimeout(r, 300); });
    const late = visible();
    t.done();
    return { early, late };
  });
  H.check('a longer task appears once it outlasts the show delay', !slow.early && slow.late, JSON.stringify(slow));
  await page.waitForTimeout(800);

  // --- a panel rebuilt mid-task picks the task up ---------------------------
  const rebuilt = await page.evaluate(async () => {
    const { startFieldTask } = await import('./state/fieldTasks.js');
    const { removeFieldPanel, addFieldPanel } = await import('./ui/FieldPanel.js');
    const t = startFieldTask('Survives rebuild', { immediate: true });
    removeFieldPanel();
    const gone = !document.querySelector('#fieldCatalogMount .field-progress');
    addFieldPanel();
    await new Promise((r) => { setTimeout(r, 300); });
    const labels = [...document.querySelectorAll('#fieldCatalogMount .field-progress-label')]
      .map((el) => el.textContent);
    t.done();
    return { gone, labels };
  });
  H.check('a rebuilt field panel shows a task that is still running',
    rebuilt.gone && rebuilt.labels.includes('Survives rebuild'), JSON.stringify(rebuilt));
  await page.waitForTimeout(800);

  // --- producers --------------------------------------------------------------
  // SCF-NCI on the blob.
  await page.evaluate(() => document.getElementById('nciScfBtn').click());
  await H.waitFor(page, () => /Added 2 fields/.test(document.getElementById('nciStatus')?.textContent || ''),
    { timeout: 60000, interval: 200 });
  // "Combine fields" on a large grid (two 104³ fields added directly).
  await page.evaluate(async () => {
    const { fieldBrowser } = await import('./ui/FieldPanel.js');
    const { Field, computeFieldStats } = await import('./model/index.js');
    const n = 104;
    const make = (label, scale) => {
      const values = new Float32Array(n * n * n);
      for (let i = 0; i < values.length; i++) values[i] = scale * Math.sin(i * 0.001);
      return new Field({
        nx: n, ny: n, nz: n, voxel: [[0.1, 0, 0], [0, 0.1, 0], [0, 0, 0.1]], values, label,
        ...computeFieldStats(values),
      });
    };
    fieldBrowser.catalog.addDerivedField(make('Big A', 1));
    fieldBrowser.catalog.addDerivedField(make('Big B', 2));
  });
  await page.waitForTimeout(500);
  await page.evaluate(() => {
    const inputs = [...document.querySelectorAll('#fieldCatalogMount .field-catalog-derived input[type="number"]')];
    // Weights follow the loaded-field order; the two big fields are the last two.
    for (const input of inputs) input.value = '0';
    inputs[inputs.length - 1].value = '1';
    inputs[inputs.length - 2].value = '1';
    document.querySelector('#fieldCatalogMount .field-catalog-derived-build').click();
  });
  await page.waitForTimeout(3000);

  const log = await page.evaluate(() => ({ tasks: [...new Set(window.__taskLog)], rows: window.__rowsSeen }));
  H.check('SCF-NCI reports a task that finishes',
    log.tasks.some((t) => /^SCF-NCI of ".*"\|done/.test(t)), JSON.stringify(log.tasks));
  H.check('combining large fields reports an immediate task that finishes',
    log.tasks.some((t) => t === 'Combining 2 fields|done|true'), JSON.stringify(log.tasks));
  H.check('the combine bar was on screen while it ran',
    log.rows.some((r) => r.startsWith('Combining 2 fields')), JSON.stringify(log.rows));

  // Loading another cube while the panel is open.
  await page.evaluate(async (text) => {
    const cv = await import('./core/crystal-viewer.js');
    await cv.loadStructure(text, 'second.cube');
  }, blobCube(40, 'Second'));
  await page.waitForTimeout(1000);
  const log2 = await page.evaluate(() => [...new Set(window.__taskLog)]);
  H.check('loading a cube file reports a task that finishes',
    log2.some((t) => t.startsWith('Loading second.cube|done')), JSON.stringify(log2));

  const end = await readBar(page);
  H.check('everything finished: the block is hidden again', end.hidden, JSON.stringify(end));
  H.check('no console errors', errors.length === 0, errors.join('\n'));

  await H.finish(browser);
})().catch(H.crash);
