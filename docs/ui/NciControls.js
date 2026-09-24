import { createNciFields } from '../model/NciField.js';
import { createInfoButton } from './InfoPanel.js';

/**
 * The "NCI analysis" block under the field list (ui/FieldPanel.js).
 *
 * Two buttons, one per way of getting the density the NCI plot is computed
 * from — finite differences on the selected field (SCF) or a sum of free-atom
 * densities on its grid (promolecular) — plus an "i" note explaining both.
 * The work itself runs in the worker pool (workers/nciTasks.js); this module
 * is only the button state, the status line and the hand-off of the two new
 * fields to the catalog.
 *
 * It lives outside FieldCatalogWidget's "Combine fields" block on purpose: that
 * block only appears with two or more loaded fields, and a Gaussian cube — the
 * file NCI is most often run on — has one.
 */

/**
 * A computation outlives the panel that started it: switching structures or
 * collapsing the panel rebuilds the DOM, and the rebuilt block must still show
 * that one is running (and must not let a second one start on top of it). So
 * the run state is module-level and every live block renders from it.
 */
const runState = {
  busy: false,
  /** @type {string} */
  status: '',
  /** @type {'info' | 'warning'} */
  tone: 'info',
  /** The catalog the status is about, so another file's panel does not show it. */
  catalog: null,
};

/** @type {Set<() => void>} render callbacks of every mounted block */
const liveBlocks = new Set();

function publish(patch) {
  Object.assign(runState, patch);
  for (const render of liveBlocks) {
    try { render(); } catch (error) { console.error('NCI controls render failed', error); }
  }
}

/** 123456 → "123 456", with a no-break space so it never wraps mid-number. */
function groupDigits(n) {
  return String(Math.round(n)).replace(/\B(?=(\d{3})+(?!\d))/g, '\u00a0');
}

/** The status line after a successful run. */
function describeResult(info, seconds) {
  const share = info.gridPoints > 0 ? (100 * info.included) / info.gridPoints : 0;
  const shareText = share > 0 && share < 0.1 ? '<0.1' : share.toFixed(1);
  const parts = [
    `Added 2 fields; ${groupDigits(info.included)} grid points in the NCI window `
    + `(${shareText}% of the grid), ${seconds.toFixed(1)} s.`,
  ];
  let warning = false;
  if (info.included === 0) {
    warning = true;
    parts.push(info.kind === 'scf'
      ? 'No point fell inside the density window, so the s = 0.5 surface is empty — '
        + 'is the selected field a total electron density?'
      : 'No point fell inside the density window, so the s = 0.5 surface is empty.');
  }
  if (info.unitWarning) {
    warning = true;
    parts.push(info.unitWarning);
  }
  if (info.clampedElements?.length) {
    warning = true;
    parts.push(`${info.clampedElements.join(', ')} approximated by the Ar free-atom density.`);
  }
  if (info.ignoredElements?.length) {
    warning = true;
    parts.push(`Ignored sites that are not elements: ${info.ignoredElements.join(', ')}.`);
  }
  return { text: parts.join(' '), tone: warning ? 'warning' : 'info' };
}

/**
 * Build the NCI block into `container`.
 *
 * @param {object} options
 * @param {HTMLElement} options.container
 * @param {() => import('../model/Field.js').Field | null} options.getSelectedField
 * @param {() => any} options.getStructure
 * @param {() => import('../model/FieldCatalog.js').FieldCatalog | null} options.getCatalog
 * @param {(field: import('../model/Field.js').Field,
 *   extra: {colourField: import('../model/Field.js').Field, info: any}) => void} options.onFieldsCreated
 *   called with the new s field once both fields are in the catalog
 * @param {(error: Error) => void} [options.onError]
 * @returns {{refresh: () => void, destroy: () => void}}
 */
export function createNciControls({
  container, getSelectedField, getStructure, getCatalog, onFieldsCreated, onError,
}) {
  container.innerHTML = '';
  const root = document.createElement('div');
  root.className = 'nci-controls';

  const heading = document.createElement('h4');
  heading.className = 'nci-controls-title';
  heading.textContent = 'Non-covalent interactions (NCI)';
  heading.appendChild(createInfoButton('./data/nciInfo.md', 'About the NCI analysis'));
  root.appendChild(heading);

  const row = document.createElement('div');
  row.className = 'nci-controls-row';

  const scfButton = document.createElement('button');
  scfButton.type = 'button';
  scfButton.id = 'nciScfBtn';
  scfButton.className = 'file-action-btn planes-calc-btn';
  scfButton.textContent = 'Create SCF-NCI field';
  row.appendChild(scfButton);

  const promoButton = document.createElement('button');
  promoButton.type = 'button';
  promoButton.id = 'nciPromolecularBtn';
  promoButton.className = 'file-action-btn planes-calc-btn';
  promoButton.textContent = 'Create Promolecular NCI field';
  row.appendChild(promoButton);

  root.appendChild(row);

  const status = document.createElement('p');
  status.id = 'nciStatus';
  status.className = 'nci-controls-status';
  status.setAttribute('role', 'status');
  status.setAttribute('aria-live', 'polite');
  root.appendChild(status);

  container.appendChild(root);

  /** The grid the promolecular density is evaluated on. */
  function promolecularGrid() {
    return getSelectedField() ?? getCatalog()?.loadedFields()?.[0] ?? null;
  }

  function render() {
    const field = getSelectedField();
    const structure = getStructure();
    const atomCount = structure?.atoms?.length ?? 0;

    const scfReady = Boolean(field?.values && field.voxel);
    scfButton.disabled = runState.busy || !scfReady;
    scfButton.title = scfReady
      ? `Compute s and sign(λ₂)ρ by finite differences on "${field.label}" (must be a total electron density)`
      : 'Select a density field in the list above first';

    const promoReady = atomCount > 0 && Boolean(promolecularGrid()?.voxel);
    promoButton.disabled = runState.busy || !promoReady;
    promoButton.title = atomCount === 0
      ? 'The structure has no atoms to build a promolecular density from'
      : promoReady
        ? 'Compute s and sign(λ₂)ρ from a sum of free-atom densities, on the selected field\'s grid'
        : 'Load a field first: the promolecular density is evaluated on its grid';

    root.setAttribute('aria-busy', String(runState.busy));
    status.textContent = runState.catalog === getCatalog() ? runState.status : '';
    status.classList.toggle('nci-controls-status-warning', runState.tone === 'warning');
  }

  async function run(kind) {
    if (runState.busy) return;
    const catalog = getCatalog();
    const structure = getStructure();
    const field = kind === 'scf' ? getSelectedField() : promolecularGrid();
    if (!catalog || !field) {
      reportError(new Error('Select a field first: NCI is computed on its grid.'));
      return;
    }

    publish({
      busy: true,
      catalog,
      status: kind === 'scf'
        ? `Computing SCF-NCI on "${field.label}"…`
        : 'Computing promolecular NCI…',
      tone: 'info',
    });
    const started = performance.now();
    try {
      const { sField, colourField, info } = await createNciFields(field, {
        kind, structure, catalog,
      });
      // The colour field first, so the s field it colours is the last entry of
      // the Derived group and the one the list scrolls to.
      catalog.addDerivedField(colourField);
      catalog.addDerivedField(sField);
      const { text, tone } = describeResult(info, (performance.now() - started) / 1000);
      publish({ busy: false, status: text, tone });
      onFieldsCreated(sField, { colourField, info });
    } catch (error) {
      publish({ busy: false, status: '', tone: 'info' });
      reportError(error instanceof Error ? error : new Error(String(error)));
    }
  }

  function reportError(error) {
    console.error('NCI analysis', error);
    if (onError) onError(error);
  }

  scfButton.addEventListener('click', () => { run('scf'); });
  promoButton.addEventListener('click', () => { run('promolecular'); });

  // Loading a band or adding a field changes what the buttons can do.
  const catalog = getCatalog();
  const unsubscribe = catalog?.subscribe ? catalog.subscribe(render) : () => {};
  liveBlocks.add(render);
  render();

  return {
    refresh: render,
    destroy() {
      unsubscribe();
      liveBlocks.delete(render);
      container.innerHTML = '';
    },
  };
}
