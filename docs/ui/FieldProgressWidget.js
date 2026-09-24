import {
  fieldTasks,
  subscribeFieldTasks,
  FIELD_TASK_SHOW_DELAY_MS,
} from '../state/fieldTasks.js';

/**
 * The progress bar under the field list: one row per running field task
 * (state/fieldTasks.js) — a label, the current stage and a bar.
 *
 * A row appears once its task has run for FIELD_TASK_SHOW_DELAY_MS (or at once
 * for an `immediate` task), creeps through the task's current stage while it
 * runs, fills and fades when it finishes, and is removed with the task. With no
 * visible task the whole block is hidden, so it takes no space.
 *
 * The widget only reads the registry: tasks come from wherever the work is
 * (FieldCatalog's on-demand loads, the NCI buttons, "Combine fields", file
 * loading), and a panel rebuilt mid-task shows the task again.
 */

/** Bar refresh interval while a task is visible. The progress is qualitative;
 *  ten updates a second keeps it moving without costing anything. */
const TICK_MS = 100;

/**
 * @param {HTMLElement} container
 * @returns {{destroy: () => void}}
 */
export function createFieldProgressWidget(container) {
  const root = document.createElement('div');
  root.className = 'field-progress';
  root.setAttribute('role', 'status');
  root.setAttribute('aria-live', 'polite');
  root.hidden = true;
  container.appendChild(root);

  /** @type {Map<number, {row: HTMLElement, text: HTMLElement, bar: HTMLProgressElement}>} */
  const rows = new Map();
  let tick = null;
  let delayTimer = null;

  function render() {
    const now = performance.now();
    const live = fieldTasks();
    const liveIds = new Set(live.map((t) => t.id));

    for (const [id, entry] of rows) {
      if (!liveIds.has(id)) {
        entry.row.remove();
        rows.delete(id);
      }
    }

    let anyVisible = false;
    let anyRunning = false;
    let pendingDelay = Infinity;
    for (const task of live) {
      if (!task.isVisible(now)) {
        if (task.state === 'running') {
          pendingDelay = Math.min(pendingDelay,
            task.startedAt + FIELD_TASK_SHOW_DELAY_MS - now);
        }
        continue;
      }
      anyVisible = true;
      if (task.state === 'running') anyRunning = true;
      let entry = rows.get(task.id);
      if (!entry) {
        entry = buildRow();
        rows.set(task.id, entry);
        root.appendChild(entry.row);
      }
      entry.text.textContent = task.stage ? `${task.label} — ${task.stage}` : task.label;
      entry.bar.value = task.displayFraction(now);
      entry.row.classList.toggle('field-progress-row-done', task.state === 'done');
      entry.row.classList.toggle('field-progress-row-failed', task.state === 'failed');
    }
    root.hidden = !anyVisible;

    if (anyRunning && tick === null) {
      tick = setInterval(render, TICK_MS);
    } else if (!anyRunning && tick !== null) {
      clearInterval(tick);
      tick = null;
    }

    // A task still inside its show delay: come back when the delay is over
    // (a task that finishes first never shows).
    if (delayTimer !== null) {
      clearTimeout(delayTimer);
      delayTimer = null;
    }
    if (Number.isFinite(pendingDelay)) {
      delayTimer = setTimeout(() => {
        delayTimer = null;
        render();
      }, Math.max(0, pendingDelay) + 5);
    }
  }

  function buildRow() {
    const row = document.createElement('div');
    row.className = 'field-progress-row';
    const text = document.createElement('span');
    text.className = 'field-progress-label';
    const bar = document.createElement('progress');
    bar.className = 'field-progress-bar';
    bar.max = 1;
    row.append(text, bar);
    return { row, text, bar };
  }

  const unsubscribe = subscribeFieldTasks(render);
  render();

  return {
    destroy() {
      unsubscribe();
      if (tick !== null) clearInterval(tick);
      if (delayTimer !== null) clearTimeout(delayTimer);
      root.remove();
      rows.clear();
    },
  };
}
