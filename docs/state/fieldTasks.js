/**
 * Long-running work on volumetric fields — reading a large field file,
 * expanding a WAVECAR band, computing derived fields (NCI, combinations) —
 * reported for the progress bar under the field list
 * (ui/FieldProgressWidget.js).
 *
 * No DOM here, so model code (FieldCatalog's on-demand loads) can report as
 * well as UI code, and a field panel rebuilt mid-task picks the task up again.
 *
 * Progress is qualitative. Most of this work is one WASM or worker call with no
 * intermediate reports, so a task is described as a sequence of stages, each
 * with the fraction range it covers; the widget creeps the bar through the
 * current range while the stage runs. A task that does get real numbers
 * reports them with setProgress().
 *
 * Nothing here is persisted.
 */

/** How long a task must run before the widget shows it, so quick jobs do not
 *  flash a bar. Tasks started with `immediate` skip this. */
export const FIELD_TASK_SHOW_DELAY_MS = 250;

/** How long a finished task stays (full, fading) before it is removed. */
export const FIELD_TASK_LINGER_MS = 600;

/** Grids at least this large count as "large": synchronous work on them
 *  (combining, marching cubes) blocks long enough to be worth a bar. */
export const LARGE_FIELD_POINTS = 1_000_000;

/** @typedef {'running' | 'done' | 'failed'} FieldTaskState */

/** @type {Map<number, FieldTask>} */
const tasks = new Map();
/** @type {Set<() => void>} */
const listeners = new Set();
let nextId = 1;

function notify() {
  for (const listener of [...listeners]) {
    try {
      listener();
    } catch (error) {
      console.error('field task listener', error);
    }
  }
}

export class FieldTask {
  /**
   * @param {string} label
   * @param {{immediate?: boolean}} [options]
   */
  constructor(label, { immediate = false } = {}) {
    this.id = nextId++;
    this.label = label;
    /** Show without waiting FIELD_TASK_SHOW_DELAY_MS — for work that is about
     *  to block the main thread, where the delay timer could not fire. */
    this.immediate = immediate;
    this.startedAt = now();
    /** @type {FieldTaskState} */
    this.state = 'running';
    /** Current stage text, e.g. 'Parsing'. */
    this.stage = '';
    /** Fraction range of the current stage, which the bar creeps through. */
    this.stageFrom = 0;
    this.stageTo = 0.9;
    this.stageStartedAt = this.startedAt;
    /** A reported fraction (setProgress), overriding the creep. @type {number | null} */
    this.fraction = null;
    /** @type {number | null} */
    this.finishedAt = null;
  }

  /**
   * Enter a stage covering [from, to] of the bar.
   * @param {string} text
   * @param {number} [from] defaults to where the bar is now
   * @param {number} [to]
   */
  setStage(text, from = this.displayFraction(), to = Math.max(from, 0.95)) {
    if (this.state !== 'running') return this;
    this.stage = text;
    this.stageFrom = clamp01(from);
    this.stageTo = Math.max(this.stageFrom, clamp01(to));
    this.stageStartedAt = now();
    this.fraction = null;
    notify();
    return this;
  }

  /** Report a real fraction in [0, 1]. */
  setProgress(fraction) {
    if (this.state !== 'running' || !Number.isFinite(fraction)) return this;
    this.fraction = clamp01(fraction);
    notify();
    return this;
  }

  done() { return this._finish('done'); }

  fail() { return this._finish('failed'); }

  _finish(state) {
    if (this.state !== 'running') return this;
    this.state = state;
    this.finishedAt = now();
    notify();
    setTimeout(() => {
      tasks.delete(this.id);
      notify();
    }, FIELD_TASK_LINGER_MS);
    return this;
  }

  /**
   * Where the bar should be at time `t`: a reported fraction if there is one;
   * otherwise the current stage's range, approached asymptotically (half-way
   * after about 1.4 s) so the bar keeps moving without ever claiming the stage
   * is finished. Full once done.
   * @param {number} [t]
   */
  displayFraction(t = now()) {
    if (this.state === 'done') return 1;
    if (this.fraction !== null) return this.fraction;
    const elapsed = Math.max(0, t - this.stageStartedAt) / 1000;
    const eased = 1 - Math.exp(-elapsed / 2);
    return this.stageFrom + (this.stageTo - this.stageFrom) * eased;
  }

  /** Whether the widget should show this task at time `t`. */
  isVisible(t = now()) {
    if (this.immediate) return true;
    const end = this.finishedAt ?? t;
    return end - this.startedAt >= FIELD_TASK_SHOW_DELAY_MS;
  }
}

/**
 * Register a task. The caller must end it with done() or fail(); prefer
 * runFieldTask(), which does that.
 * @param {string} label
 * @param {{immediate?: boolean}} [options]
 */
export function startFieldTask(label, options) {
  const task = new FieldTask(label, options);
  tasks.set(task.id, task);
  notify();
  return task;
}

/**
 * Run `work` as a task: done() when it resolves, fail() when it throws (the
 * error is re-thrown).
 * @template T
 * @param {string} label
 * @param {(task: FieldTask) => Promise<T> | T} work
 * @param {{immediate?: boolean}} [options]
 * @returns {Promise<T>}
 */
export async function runFieldTask(label, work, options) {
  const task = startFieldTask(label, options);
  try {
    const result = await work(task);
    task.done();
    return result;
  } catch (error) {
    task.fail();
    throw error;
  }
}

/** Live tasks, oldest first (including finished ones still lingering). */
export function fieldTasks() {
  return [...tasks.values()];
}

/**
 * @param {() => void} listener called on every task change
 * @returns {() => void} unsubscribe
 */
export function subscribeFieldTasks(listener) {
  listeners.add(listener);
  return () => { listeners.delete(listener); };
}

/** Whether synchronous work on this field is slow enough to show a bar for. */
export function isLargeField(field) {
  const n = (field?.nx ?? 0) * (field?.ny ?? 0) * (field?.nz ?? 0);
  return n >= LARGE_FIELD_POINTS;
}

/**
 * Resolve after the browser has painted, so a bar shown just now is on screen
 * before synchronous work blocks the main thread. Two animation frames: the
 * first callback runs before the paint, the second after it. Falls back to a
 * timeout where there are no animation frames (node, a hidden tab).
 * @returns {Promise<void>}
 */
export function nextPaint() {
  return new Promise((resolve) => {
    let settled = false;
    const finish = () => {
      if (!settled) {
        settled = true;
        resolve();
      }
    };
    setTimeout(finish, 50);
    if (typeof requestAnimationFrame === 'function') {
      requestAnimationFrame(() => requestAnimationFrame(finish));
    }
  });
}

function now() {
  return typeof performance !== 'undefined' ? performance.now() : Date.now();
}

function clamp01(x) {
  return Math.min(1, Math.max(0, Number(x) || 0));
}
