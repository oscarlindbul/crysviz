// "Colour surface by" control of the Volumetric Field panel (ui/FieldPanel.js):
// colours the selected field's isosurface per vertex by another loaded field
// on the same grid, through a colormap over [min, max] — the NCI plot's
// s = 0.5 surface coloured by sign(λ₂)ρ is the motivating case, but any pair
// of same-grid fields works.
//
// The state lives on the Field itself (Field.colorBy, see model/Field.js), so
// a colour-by set programmatically (the NCI buttons preset one on the field
// they create) is simply read back by sync(). Nothing here is persisted to
// browser storage: derived fields are not persisted either.

import { groups } from '../state/store.js';
import { updateField, requestRender } from '../render/index.js';
import { canColorFieldBy, refreshIsosurfaceColors, resolveColorBy } from '../model/Isosurface.js';
import { createColorBar, COLORMAP_CHOICES } from './ColorBarWidget.js';
import { registerColorBarSource } from './ColorBarRegistry.js';

const BAR_FLOATING_ID = 'fieldColorByBarFloating';
const DEFAULT_COLORMAP = 'bgyor';

/** The live colour bar (null while colour-by is off), for image export. */
let barInstance = null;
registerColorBarSource('isocolor', 'Isosurface colour', () => barInstance);

/** Round for display without shrinking a range to nothing. */
function tidy(value) {
  return Number(Number(value).toPrecision(4));
}

/**
 * The starting range for colouring by `colorField`: symmetric ±max|v| for a
 * field that takes both signs (so 0 sits mid-map, as a diverging map like
 * BGYOR expects), else the field's own [min, max].
 *
 * @param {import('../model/Field.js').Field} colorField
 * @returns {{min: number, max: number}}
 */
export function defaultColorByRange(colorField) {
  let lo = Number(colorField?.minValue);
  let hi = Number(colorField?.maxValue);
  if (!Number.isFinite(lo) || !Number.isFinite(hi)) {
    lo = Infinity;
    hi = -Infinity;
    const values = colorField?.values ?? [];
    for (let i = 0; i < values.length; i++) {
      const v = values[i];
      if (v < lo) lo = v;
      if (v > hi) hi = v;
    }
    if (!Number.isFinite(lo) || !Number.isFinite(hi)) return { min: -1, max: 1 };
  }
  if (lo < 0 && hi > 0) {
    const a = tidy(Math.max(-lo, hi));
    return { min: -a, max: a };
  }
  if (!(hi > lo)) return { min: tidy(lo) - 1, max: tidy(lo) + 1 };
  return { min: tidy(lo), max: tidy(hi) };
}

/**
 * Build the control into `container` (its children are replaced).
 *
 * @param {HTMLElement} container
 * @param {{ getField: () => (import('../model/Field.js').Field | null),
 *   getCandidates: () => import('../model/Field.js').Field[],
 *   catalog?: { subscribe: (listener: () => void) => () => void } | null }} options
 *   getField: the selected field (whose isosurface is coloured);
 *   getCandidates: every loaded field (filtered here to the same grid);
 *   catalog: re-sync when fields are added or evicted.
 * @returns {{ sync: () => void, destroy: () => void }}
 */
export function createFieldColorByControl(container, { getField, getCandidates, catalog = null }) {
  container.innerHTML = '';
  container.classList.add('field-colorby');

  const label = document.createElement('label');
  label.htmlFor = 'fieldColorBySelect';
  label.textContent = 'Colour surface by:';
  const fieldSelect = document.createElement('select');
  fieldSelect.id = 'fieldColorBySelect';
  fieldSelect.className = 'planes-select planes-full-width';
  fieldSelect.title = 'Colour the isosurface by the values of another field on the same grid';

  const options = document.createElement('div');
  options.id = 'fieldColorByOptions';
  options.className = 'field-colorby-options';
  options.hidden = true;

  const cmapLabel = document.createElement('label');
  cmapLabel.htmlFor = 'fieldColorByColormap';
  cmapLabel.textContent = 'Colormap:';
  const cmapSelect = document.createElement('select');
  cmapSelect.id = 'fieldColorByColormap';
  cmapSelect.className = 'planes-select planes-full-width';
  for (const [value, text] of COLORMAP_CHOICES) {
    const opt = document.createElement('option');
    opt.value = value;
    opt.textContent = text;
    cmapSelect.appendChild(opt);
  }
  cmapSelect.value = DEFAULT_COLORMAP;

  // The bar carries the Min/Max inputs (its onLimitsCommit) and the legend.
  const barMount = document.createElement('div');
  barMount.id = 'fieldColorByBarMount';
  barMount.className = 'field-colorby-bar';

  options.append(cmapLabel, cmapSelect, barMount);
  container.append(label, fieldSelect, options);

  /** Fields offered in the dropdown, by option index. */
  let candidates = [];
  /** The colour field the current bar was built for (its legend). */
  let barField = null;

  function removeBar() {
    barInstance?.remove();
    barInstance = null;
    barField = null;
  }

  /** Push the field's colour-by into the scene: remap when possible, else
   *  rebuild (a different colour field needs a new marching-cubes pass). */
  function applyToScene(field, { rebuild = false } = {}) {
    if (!field || groups.activeField !== field || !groups.isosurfaceGroup) return;
    if (rebuild || !refreshIsosurfaceColors()) updateField();
    requestRender();
  }

  function setColorBy(field, colorBy, opts) {
    field.colorBy = colorBy;
    applyToScene(field, opts);
    sync();
  }

  function syncBar(colorBy) {
    if (!colorBy) {
      removeBar();
      return;
    }
    if (barInstance && barField === colorBy.field) {
      barInstance.update(colorBy.colormap);
      barInstance.setRange(colorBy.min, colorBy.max);
      return;
    }
    removeBar();
    barField = colorBy.field;
    const fallback = defaultColorByRange(colorBy.field);
    barInstance = createColorBar(barMount, colorBy.colormap, colorBy.min, colorBy.max, {
      floatingId: BAR_FLOATING_ID,
      fallbackMin: fallback.min,
      fallbackMax: fallback.max,
      legend: colorBy.field.label || 'Colour field',
      onLimitsCommit: (min, max) => {
        const field = getField();
        if (!field?.colorBy) return;
        setColorBy(field, { ...field.colorBy, min, max });
      },
      onAutoRange: () => {
        const field = getField();
        if (!field?.colorBy) return;
        setColorBy(field, { ...field.colorBy, ...defaultColorByRange(field.colorBy.field) });
      },
    });
  }

  function sync() {
    const field = getField();
    // A colour-by that cannot be honoured (grid mismatch, no values) reads
    // as none here, exactly as the isosurface treats it.
    const colorBy = field ? resolveColorBy(field) : null;
    candidates = field
      ? getCandidates().filter((f) => f !== field && canColorFieldBy(field, f))
      : [];
    // Keep a programmatic choice visible even when the list would not offer
    // it (the field itself, or a derived field since removed from the list).
    if (colorBy && !candidates.includes(colorBy.field)) candidates.push(colorBy.field);

    fieldSelect.replaceChildren();
    const none = document.createElement('option');
    none.value = '';
    none.textContent = 'None (flat colour)';
    fieldSelect.appendChild(none);
    candidates.forEach((f, i) => {
      const opt = document.createElement('option');
      opt.value = String(i);
      opt.textContent = f === field ? `${f.label || 'This field'} (itself)` : (f.label || `Field ${i + 1}`);
      fieldSelect.appendChild(opt);
    });
    fieldSelect.value = colorBy ? String(candidates.indexOf(colorBy.field)) : '';
    fieldSelect.disabled = !field || candidates.length === 0;

    options.hidden = !colorBy;
    if (colorBy) cmapSelect.value = colorBy.colormap;
    syncBar(colorBy);
  }

  fieldSelect.addEventListener('change', () => {
    const field = getField();
    if (!field) return;
    const colorField = fieldSelect.value === '' ? null : candidates[Number(fieldSelect.value)];
    if (!colorField) {
      setColorBy(field, null);
      return;
    }
    const colormap = field.colorBy?.colormap || cmapSelect.value || DEFAULT_COLORMAP;
    setColorBy(field, { field: colorField, colormap, ...defaultColorByRange(colorField) }, { rebuild: true });
  });

  cmapSelect.addEventListener('change', () => {
    const field = getField();
    if (!field?.colorBy) return;
    setColorBy(field, { ...field.colorBy, colormap: cmapSelect.value });
  });

  const unsubscribe = catalog?.subscribe(() => sync()) ?? null;
  sync();

  return {
    sync,
    destroy() {
      unsubscribe?.();
      removeBar();
      container.innerHTML = '';
    },
  };
}
