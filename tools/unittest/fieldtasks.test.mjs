// The field-task registry behind the progress bar under the field list
// (docs/state/fieldTasks.js). No DOM: the widget is covered by
// tools/browsertest/tests/fieldprogress.test.js.
import { test } from 'node:test';
import assert from 'node:assert/strict';

import {
  startFieldTask,
  runFieldTask,
  fieldTasks,
  subscribeFieldTasks,
  isLargeField,
  FIELD_TASK_SHOW_DELAY_MS,
  FIELD_TASK_LINGER_MS,
  LARGE_FIELD_POINTS,
} from '../../docs/state/fieldTasks.js';

const sleep = (ms) => new Promise((resolve) => { setTimeout(resolve, ms); });

test('a task is listed while running, full when done, and removed after lingering', async () => {
  const task = startFieldTask('Loading x');
  assert.ok(fieldTasks().includes(task));
  assert.equal(task.state, 'running');
  task.done();
  assert.equal(task.state, 'done');
  assert.equal(task.displayFraction(), 1);
  assert.ok(fieldTasks().includes(task), 'still listed while it fades');
  await sleep(FIELD_TASK_LINGER_MS + 50);
  assert.ok(!fieldTasks().includes(task));
});

test('the bar creeps through a stage without leaving its range', () => {
  const task = startFieldTask('Stage test');
  task.setStage('computing', 0.2, 0.6);
  const t0 = task.stageStartedAt;
  const early = task.displayFraction(t0 + 10);
  const later = task.displayFraction(t0 + 3000);
  const much = task.displayFraction(t0 + 1e6);
  assert.ok(early >= 0.2 && early < 0.25, `early ${early}`);
  assert.ok(later > early && later < 0.6, `later ${later}`);
  assert.ok(much <= 0.6 && much > 0.59, `much later ${much}`);
  // A later stage defaults to starting where the bar is.
  task.setStage('meshing', undefined, 0.9);
  assert.ok(task.stageFrom >= 0.2 && task.stageTo === 0.9);
  // Reported progress overrides the creep until the next stage.
  task.setProgress(0.75);
  assert.equal(task.displayFraction(t0 + 1e6), 0.75);
  task.done();
});

test('only tasks that outlast the show delay (or are immediate) are visible', () => {
  const quick = startFieldTask('quick');
  const t = quick.startedAt;
  assert.equal(quick.isVisible(t + FIELD_TASK_SHOW_DELAY_MS - 1), false);
  assert.equal(quick.isVisible(t + FIELD_TASK_SHOW_DELAY_MS + 1), true);
  quick.finishedAt = t + 10;
  assert.equal(quick.isVisible(t + 1000), false, 'finished inside the delay: never shown');
  quick.finishedAt = null;
  quick.done();

  const now = startFieldTask('blocking', { immediate: true });
  assert.equal(now.isVisible(now.startedAt), true);
  now.done();
});

test('runFieldTask ends the task either way and rethrows', async () => {
  const seen = [];
  const unsubscribe = subscribeFieldTasks(() => {
    for (const t of fieldTasks()) seen.push(`${t.label}:${t.state}`);
  });
  const value = await runFieldTask('ok', async () => 42);
  assert.equal(value, 42);
  await assert.rejects(runFieldTask('bad', async () => { throw new Error('boom'); }), /boom/);
  unsubscribe();
  assert.ok(seen.includes('ok:done'));
  assert.ok(seen.includes('bad:failed'));
});

test('finished tasks ignore further updates', () => {
  const task = startFieldTask('closed');
  task.fail();
  task.setStage('late', 0, 1);
  task.setProgress(0.5);
  task.done();
  assert.equal(task.state, 'failed');
  assert.equal(task.stage, '');
  assert.equal(task.fraction, null);
});

test('isLargeField uses the grid size', () => {
  assert.equal(isLargeField({ nx: 100, ny: 100, nz: 99 }), false);
  assert.equal(isLargeField({ nx: 100, ny: 100, nz: 100 }), 100 * 100 * 100 >= LARGE_FIELD_POINTS);
  assert.equal(isLargeField({ nx: 200, ny: 200, nz: 200 }), true);
  assert.equal(isLargeField(null), false);
});

test('an on-demand catalog load (a WAVECAR band) is reported as a task', async () => {
  const { FieldCatalogNode } = await import('../../docs/model/FieldCatalog.js');
  const node = new FieldCatalogNode({
    id: 'wf:1:3:12',
    label: 'Band 12',
    meta: { spin: 1, kpt: 3, band: 12 },
    load: async () => ({ label: 'loaded' }),
  });
  const seen = new Set();
  const unsubscribe = subscribeFieldTasks(() => {
    for (const t of fieldTasks()) seen.add(`${t.label}|${t.state}`);
  });
  const field = await node.ensureLoaded();
  unsubscribe();
  assert.equal(field.label, 'loaded');
  assert.ok(seen.has('Loading k-point 3, band 12|running'), [...seen].join(', '));
  assert.ok(seen.has('Loading k-point 3, band 12|done'), [...seen].join(', '));
});
