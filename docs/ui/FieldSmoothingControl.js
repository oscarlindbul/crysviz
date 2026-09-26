// "Smoothing" control of the Volumetric Field panel (ui/FieldPanel.js): one
// drop-down picking the active isosurface smoothing method, and under it the
// parameter controls of that method, rebuilt from its spec
// (SMOOTHING_METHODS in model/Isosurface.js) whenever the method changes.
//
// The settings are GLOBAL module state in model/Isosurface.js, like the
// isosurface colours/opacity, and each method keeps its own values there, so
// switching back to a method restores what it was set to. This control only
// edits them and tells the panel (onChange) to re-mesh; it keeps no copy.

import {
  SMOOTHING_METHODS,
  getIsosurfaceSmoothingSettings,
  setIsosurfaceSmoothingSettings,
} from '../model/Isosurface.js';
import { createInfoButton } from './InfoPanel.js';

/** The live control, so a state restore can re-read the settings into it. */
let activeControl = null;

/** Re-read the global smoothing settings into the mounted control (if any),
 *  e.g. after ui/ShareModule.js restored them from a saved session. */
export function syncFieldSmoothingControl() {
  activeControl?.sync();
}

/** Number of decimals a readout needs for `step` (0.05 -> 2). */
function decimalsFor(step) {
  const text = String(step ?? 1);
  const dot = text.indexOf('.');
  return dot < 0 ? 0 : text.length - dot - 1;
}

/**
 * Build the control into `container` (its children are replaced).
 *
 * @param {HTMLElement} container
 * @param {{ onChange: (info: {live: boolean}) => void }} options
 *   onChange: the settings changed; `live` is true while a slider is being
 *   dragged (the panel coalesces those) and false on a committed change.
 * @returns {{ sync: () => void, destroy: () => void }}
 */
export function createFieldSmoothingControl(container, { onChange }) {
  container.innerHTML = '';
  container.classList.add('field-smoothing');

  // One row: label, method drop-down, (i). The explanations of the methods
  // and their settings live in the info panel, not in the control.
  const head = document.createElement('div');
  head.className = 'field-smoothing-head';
  const label = document.createElement('label');
  label.htmlFor = 'FieldSmoothingMethod';
  label.textContent = 'Smoothing:';
  const methodSelect = document.createElement('select');
  methodSelect.id = 'FieldSmoothingMethod';
  methodSelect.className = 'planes-select';
  for (const method of SMOOTHING_METHODS) {
    const opt = document.createElement('option');
    opt.value = method.id;
    opt.textContent = method.label;
    methodSelect.appendChild(opt);
  }
  head.append(label, methodSelect, createInfoButton('./data/fieldSmoothingInfo.md', 'About surface smoothing'));

  // label | slider | readout grid, one row per parameter
  const paramsBox = document.createElement('div');
  paramsBox.id = 'FieldSmoothingParams';
  paramsBox.className = 'field-smoothing-params';

  container.append(head, paramsBox);

  /** Write one parameter of the active method and notify. */
  function setParam(methodId, key, value, live) {
    setIsosurfaceSmoothingSettings({ params: { [methodId]: { [key]: value } } });
    onChange({ live });
  }

  /** One row per parameter of `method`, holding `values`. Element ids are
   *  FieldSmoothing-<key> (plus FieldSmoothing-<key>-value for a readout). */
  function buildParams(method, values) {
    paramsBox.replaceChildren();
    for (const spec of method.params) {
      const id = `FieldSmoothing-${spec.key}`;
      const value = values[spec.key] ?? spec.default;

      if (spec.type === 'bool') {
        const row = document.createElement('label');
        row.className = 'toggle_row toggle_container field-smoothing-wide-row';
        const sw = document.createElement('span');
        sw.className = 'toggle_switch';
        const input = document.createElement('input');
        input.type = 'checkbox';
        input.id = id;
        input.checked = Boolean(value);
        const knob = document.createElement('span');
        knob.className = 'toggle_slider';
        sw.append(input, knob);
        const text = document.createElement('span');
        text.className = 'toggle_text';
        text.textContent = ` ${spec.label}`;
        row.append(sw, text);
        input.addEventListener('change', () => setParam(method.id, spec.key, input.checked, false));
        paramsBox.appendChild(row);
        continue;
      }

      const rowLabel = document.createElement('label');
      rowLabel.htmlFor = id;
      rowLabel.textContent = `${spec.label}:`;

      if (spec.type === 'choice') {
        const select = document.createElement('select');
        select.id = id;
        select.className = 'planes-select field-smoothing-wide';
        for (const choice of spec.choices ?? []) {
          const opt = document.createElement('option');
          opt.value = choice.value;
          opt.textContent = choice.label;
          select.appendChild(opt);
        }
        select.value = String(value);
        select.addEventListener('change', () => setParam(method.id, spec.key, select.value, false));
        paramsBox.append(rowLabel, select);
        continue;
      }

      // int / float: range + numeric readout in the grid's last two columns
      const range = document.createElement('input');
      range.type = 'range';
      range.id = id;
      range.min = String(spec.min);
      range.max = String(spec.max);
      range.step = String(spec.step ?? (spec.type === 'int' ? 1 : 0.01));
      range.value = String(value);
      const readout = document.createElement('span');
      readout.id = `${id}-value`;
      const decimals = spec.type === 'int' ? 0 : decimalsFor(spec.step);
      const show = (v) => { readout.textContent = Number(v).toFixed(decimals); };
      show(value);
      // Drag: live (coalesced by the panel); release: the committed value.
      range.addEventListener('input', () => {
        show(range.value);
        setParam(method.id, spec.key, Number(range.value), true);
      });
      range.addEventListener('change', () => setParam(method.id, spec.key, Number(range.value), false));
      paramsBox.append(rowLabel, range, readout);
    }
    paramsBox.hidden = method.params.length === 0;
  }

  function sync() {
    const settings = getIsosurfaceSmoothingSettings();
    const method = SMOOTHING_METHODS.find((m) => m.id === settings.method) ?? SMOOTHING_METHODS[0];
    methodSelect.value = method.id;
    buildParams(method, settings.params[method.id] ?? {});
  }

  methodSelect.addEventListener('change', () => {
    setIsosurfaceSmoothingSettings({ method: methodSelect.value });
    sync();
    onChange({ live: false });
  });

  sync();

  const control = {
    sync,
    destroy() {
      if (activeControl === control) activeControl = null;
      container.innerHTML = '';
    },
  };
  activeControl = control;
  return control;
}
