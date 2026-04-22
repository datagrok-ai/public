import * as ui from 'datagrok-api/ui';
import {PluginUIContext} from 'molstar/lib/mol-plugin-ui/context';
import {
  Structure,
  StructureElement,
  StructureProperties as SP,
  StructureSelection,
  Queries as Q,
  StructureQuery,
  Unit,
} from 'molstar/lib/mol-model/structure';
import {MolScriptBuilder as MS} from 'molstar/lib/mol-script/language/builder';
import {Expression} from 'molstar/lib/mol-script/language/expression';
import {StructureUniqueSubsetBuilder} from 'molstar/lib/mol-model/structure/structure/util/unique-subset-builder';
import {_package} from '../../package';

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

const BINDING_SITE_COMPONENT_KEY = 'datagrok-binding-site';

/** Residue names that are never considered as the ligand seed (water, common ions, cryoprotectants). */
const HET_BLOCKLIST = new Set([
  'HOH', 'WAT', 'DOD',
  'NA', 'CL', 'MG', 'CA', 'ZN', 'K', 'FE', 'MN', 'CU', 'NI', 'CO',
  'SO4', 'PO4', 'GOL', 'EDO', 'DMS', 'PEG', 'MPD', 'ACE', 'FMT',
  'IMD', 'IOD', 'BR',
]);

// ---------------------------------------------------------------------------
// Path A helpers — single Structure, HET inside PDB / mmCIF
// ---------------------------------------------------------------------------

/**
 * Returns a MolScript Expression that selects the "best" ligand HET group
 * within the given Structure, or null if no HET group is found.
 *
 * The "best" HET is the non-polymer, non-water group with the most heavy atoms,
 * excluding residue names in HET_BLOCKLIST.
 * @param {Structure} structure the Mol* Structure to scan for a ligand.
 * @return {Expression | null} the selection Expression, or null when no viable HET was found.
 */
export function getLigandSeedExpression(structure: Structure): Expression | null {
  const candidates: Array<{compId: string; chainId: string; seqId: number; count: number}> = [];

  for (const unit of structure.units) {
    if (!Unit.isAtomic(unit)) continue;
    const {elements} = unit;
    if (elements.length === 0) continue;

    const loc = StructureElement.Location.create(structure, unit, elements[0]);

    let prevResIdx = -1;
    let curCount = 0;
    let curCompId = '';
    let curChainId = '';
    let curSeqId = 0;

    const commit = () => {
      if (curCompId && !HET_BLOCKLIST.has(curCompId) && curCount > 0)
        candidates.push({compId: curCompId, chainId: curChainId, seqId: curSeqId, count: curCount});
    };

    for (let i = 0, n = elements.length; i < n; i++) {
      StructureElement.Location.set(loc, structure, unit, elements[i]);
      const entityType = SP.entity.type(loc);
      if (entityType !== 'non-polymer') continue;

      const resIdx = unit.residueIndex[elements[i]];
      if (resIdx !== prevResIdx) {
        if (prevResIdx !== -1) commit();
        prevResIdx = resIdx;
        curCompId = SP.atom.label_comp_id(loc);
        curChainId = SP.chain.label_asym_id(loc);
        curSeqId = SP.residue.label_seq_id(loc);
        curCount = 0;
      }
      const element = SP.atom.type_symbol(loc);
      if (element !== 'H' && element !== 'D')
        curCount++;
    }
    if (prevResIdx !== -1) commit();
  }

  if (candidates.length === 0) return null;

  // Pick the HET with the most heavy atoms.
  const best = candidates.reduce((a, b) => b.count > a.count ? b : a);

  return MS.struct.generator.atomGroups({
    'chain-test': MS.core.rel.eq([MS.ammp('label_asym_id'), best.chainId]),
    'residue-test': MS.core.logic.and([
      MS.core.rel.eq([MS.ammp('label_seq_id'), best.seqId]),
      MS.core.rel.eq([MS.ammp('label_comp_id'), best.compId]),
    ]),
  });
}

// ---------------------------------------------------------------------------
// Path B helpers — cross-Structure, ligand in separate trajectory
// ---------------------------------------------------------------------------

/**
 * Returns world-space atom positions from a ligand Structure's units as
 * an interleaved Float32Array [x0,y0,z0, x1,y1,z1, ...].
 * Returns null if the ligand structure has no atoms.
 * @param {Structure} ligandStructure the ligand-only Structure to read coordinates from.
 * @return {Float32Array | null} interleaved xyz positions, or null for an empty structure.
 */
export function getLigandAtomPositions(ligandStructure: Structure): Float32Array | null {
  // Pass 1: total atom count — lets us allocate the Float32Array once
  // at the final size and write via index assignment. Cheaper than growing
  // a number[] via push and re-copying into a typed array at the end.
  let atomCount = 0;
  for (const unit of ligandStructure.units)
    atomCount += unit.elements.length;
  if (atomCount === 0) return null;

  const positions = new Float32Array(atomCount * 3);
  let w = 0;
  for (const unit of ligandStructure.units) {
    const {elements, conformation} = unit;
    for (let i = 0, n = elements.length; i < n; i++) {
      const eIdx = elements[i];
      positions[w++] = conformation.x(eIdx);
      positions[w++] = conformation.y(eIdx);
      positions[w++] = conformation.z(eIdx);
    }
  }
  return positions;
}

/**
 * Builds a binding-site component on `proteinStructureRef` using a spatial query
 * against the provided `ligandPositions` array.
 * @param {PluginUIContext} plugin the Mol* plugin context that owns the state tree.
 * @param {string} proteinStructureRef state-tree ref of the protein Structure to decorate.
 * @param {Float32Array} ligandPositions interleaved xyz world coordinates defining the query sphere(s).
 * @param {number} radius per-atom sphere radius in Å.
 * @param {boolean} wholeResidues when true, expand atom hits to whole residues.
 * @return {Promise<string[]>} the refs of the created [component, representation] cells (empty on failure).
 */
export async function buildBindingSiteComponentFromPositions(
  plugin: PluginUIContext,
  proteinStructureRef: string,
  ligandPositions: Float32Array,
  radius: number,
  wholeResidues: boolean,
): Promise<string[]> {
  const cell = plugin.state.data.selectQ((q) => q.byRef(proteinStructureRef))[0];
  const proteinStructureData = cell?.obj?.data as Structure | undefined;
  if (!proteinStructureData) {
    _package.logger.warning('buildBindingSiteComponentFromPositions: protein structure not found');
    return [];
  }

  const builder = new StructureUniqueSubsetBuilder(proteinStructureData);
  for (let i = 0; i < ligandPositions.length; i += 3) {
    proteinStructureData.lookup3d.findIntoBuilder(
      ligandPositions[i], ligandPositions[i + 1], ligandPositions[i + 2],
      radius, builder,
    );
  }

  if (builder.isEmpty) return [];

  let subStructure = builder.getStructure();

  /* eslint-disable new-cap -- Mol* static factory methods use PascalCase by convention */
  if (wholeResidues) {
    try {
      const sel = StructureSelection.Singletons(proteinStructureData, subStructure);
      // Mol*'s `wholeResidues` expects a StructureQuery `(ctx) => selection`; we
      // pass a thunk that ignores ctx and returns an already-computed selection.
      // The runtime contract is respected — only the signature doesn't match.
      const thunk = (() => sel) as unknown as StructureQuery;
      const expanded = Q.modifiers.wholeResidues(thunk);
      const expandedSel = StructureQuery.run(expanded, proteinStructureData);
      subStructure = StructureSelection.unionStructure(expandedSel);
    } catch (e) {
      // If whole-residue expansion fails, fall back to the atom-only subset.
      _package.logger.warning(`wholeResidues expansion failed: ${String(e)}`);
    }
  }

  const subSel = StructureSelection.Singletons(proteinStructureData, subStructure);
  const loci = StructureSelection.toLociWithSourceUnits(subSel);
  const bundle = StructureElement.Bundle.fromLoci(loci);
  const expression = StructureElement.Bundle.toExpression(bundle);
  /* eslint-enable new-cap */

  return buildBindingSiteComponentFromExpression(plugin, proteinStructureRef, expression);
}

/**
 * Shared final step: create named component + ball-and-stick representation
 * from an Expression on proteinStructureRef.
 * @param {PluginUIContext} plugin the Mol* plugin context that owns the state tree.
 * @param {string} proteinStructureRef state-tree ref of the protein Structure to decorate.
 * @param {Expression} expression the MolScript Expression selecting the binding-site atoms.
 * @return {Promise<string[]>} the refs of the created [component, representation] cells.
 */
export async function buildBindingSiteComponentFromExpression(
  plugin: PluginUIContext,
  proteinStructureRef: string,
  expression: Expression,
): Promise<string[]> {
  const comp = await plugin.builders.structure.tryCreateComponentFromExpression(
    proteinStructureRef,
    expression,
    BINDING_SITE_COMPONENT_KEY,
    {label: 'Binding site'},
  );
  if (!comp) return [];

  const repr = await plugin.builders.structure.representation.addRepresentation(
    comp,
    {type: 'ball-and-stick'},
  );

  return [comp.ref, repr.ref];
}

// ---------------------------------------------------------------------------
// UI — right-strip overlay
// ---------------------------------------------------------------------------

export type BindingSiteOverlayHandlers = {
  getShowBindingSite: () => boolean;
  setShowBindingSite: (v: boolean) => void;
  getBindingSiteRadius: () => number;
  setBindingSiteRadius: (v: number) => void;
  isLigandAvailable: () => boolean;
  /** True if any loaded structure contains a polymer (protein). When false,
   *  the overlay is hidden entirely — side-chain visualization isn't
   *  meaningful without a polymer. */
  hasPolymer: () => boolean;
};

/** Overlay root element with the imperative methods the viewer uses to drive
 *  layout and cleanup. The methods are attached as expandos in
 *  `createBindingSiteOverlay`; declaring them here keeps callers type-safe. */
export type BindingSiteOverlayElement = HTMLElement & {
  _bsvReposition?: () => void;
  _bsvRefreshUI?: () => void;
  _bsvCleanup?: () => void;
};

/**
 * Computes the vertical offset from `ancestor`'s padding box to `descendant`'s
 * top edge, using `offsetTop`/`offsetParent` chain. Returns local (unzoomed)
 * pixels — matching what `style.top` expects — and is therefore robust to CSS
 * `zoom` applied to any ancestor. `getBoundingClientRect` returns zoom-scaled
 * pixels, which mismatches `style.top` when zoom is non-1.
 * @param {HTMLElement} descendant the element whose top edge we want.
 * @param {HTMLElement} ancestor the element to measure relative to.
 * @return {number | null} the local-pixel offset, or null if `ancestor` isn't on the offsetParent chain.
 */
function offsetTopWithin(descendant: HTMLElement, ancestor: HTMLElement): number | null {
  let y = 0;
  let el: HTMLElement | null = descendant;
  while (el && el !== ancestor) {
    y += el.offsetTop;
    el = el.offsetParent as HTMLElement | null;
  }
  return el === ancestor ? y : null; // ancestor not on offset chain
}

/**
 * Positions the overlay inside its parent (ideally `.msp-viewport`) so its
 * top edge sits just below the Mol* icon strip when present, or at the
 * parent's top-right corner otherwise. The overlay's `right` is kept at
 * 10 px so it aligns with the strip's right edge (strip is also at right:10).
 * @param {HTMLElement} wrapper the overlay element; its `top` / `right` are updated in place.
 */
function positionOverlayRelativeToStrip(wrapper: HTMLElement): void {
  const container = wrapper.parentElement;
  if (!container) return;

  // Anchor to `.msp-viewport-controls-buttons` (the inner button stack) —
  // its bottom edge matches the last icon exactly. The outer
  // `.msp-viewport-controls` wrapper has ~4 px extra padding that would
  // double our inter-icon gap. Scoped to the overlay's ancestor chain so
  // a second Biostructure viewer's strip can't be picked up by mistake.
  let controls: HTMLElement | null = null;
  let search: HTMLElement | null = container;
  while (search && !controls) {
    const found = search.querySelector('.msp-viewport-controls-buttons') as HTMLElement | null;
    if (found && found.offsetHeight > 0 && search.contains(found)) controls = found;
    search = search.parentElement;
    if (search && (search.classList?.contains('d4-root') || search === document.body)) break;
  }

  if (controls && controls.offsetHeight > 0) {
    // Prefer offsetTop chain (local pixels, zoom-safe). Fall back to rect
    // math only if the strip isn't on the overlay's offsetParent chain.
    const stripTopInContainer = offsetTopWithin(controls, container as HTMLElement);
    const newTop = stripTopInContainer !== null ?
      stripTopInContainer + controls.offsetHeight + 4 :
      Math.round(controls.getBoundingClientRect().bottom - container.getBoundingClientRect().top + 4);
    wrapper.style.top = `${newTop}px`;
    return;
  }

  // Fallback: sensible default top that's visible.
  wrapper.style.top = '200px';
}

/**
 * Creates the binding-site overlay DOM element.
 * The element is positioned absolute to visually extend the .msp-viewport-controls
 * strip downward. The caller appends it to viewer.root, then should call
 * `wrapper._bsvReposition()` once the Mol* strip is in the DOM.
 * @param {BindingSiteOverlayHandlers} handlers state getters/setters and predicates owned by the viewer.
 * @return {BindingSiteOverlayElement} the overlay root element with
 *   expando methods `_bsvReposition` / `_bsvRefreshUI` / `_bsvCleanup`.
 */
export function createBindingSiteOverlay(handlers: BindingSiteOverlayHandlers): BindingSiteOverlayElement {
  const wrapper = document.createElement('div') as BindingSiteOverlayElement;
  wrapper.className = 'bsv-bs-overlay';

  // Mol*'s native strip buttons each sit over a `.msp-semi-transparent-background`
  // sibling (#eeece7 @ 0.5 opacity). Reuse the same native class to match
  // exactly — anything smaller/larger makes our bullseye stand out.
  const iconBg = ui.div([], 'msp-semi-transparent-background');

  // Keep raw `<button>` (not `ui.button`) so we can apply Mol*'s native
  // strip-icon classes, which drive the color states for free:
  //   msp-btn-link-toggle-off → color: #9c835f (muted, default)
  //   :hover                  → color: #ae5d04 (orange accent)
  //   msp-btn-link-toggle-on  → color: #332b1f (fully saturated on active)
  // The glyph itself comes from `ui.iconFA('bullseye')` — the FA icon
  // inherits `currentColor` from the parent button automatically.
  const iconBtn = document.createElement('button');
  iconBtn.className = 'msp-btn msp-btn-icon msp-btn-link msp-btn-link-toggle-off bsv-bs-icon-btn';
  iconBtn.type = 'button';
  iconBtn.setAttribute('aria-label', 'Binding site');
  iconBtn.setAttribute('title', 'Binding site');
  iconBtn.appendChild(ui.iconFA('bullseye'));

  const iconHolder = ui.div([iconBg, iconBtn], 'bsv-bs-icon-holder');
  wrapper.appendChild(iconHolder);

  // Popover attached to wrapper.parentElement (viewer root). Width / top /
  // left / right are set at runtime by _bsvReposition based on the viewport
  // clip box — the rest lives in CSS (`.bsv-bs-popover`).
  const popover = ui.div([], 'msp-viewport-controls-panel bsv-bs-popover bsv-bs-popover-hidden');
  popover.setAttribute('role', 'dialog');
  popover.setAttribute('aria-label', 'Binding site options');

  // FA bullseye icon — size + colour driven by CSS (`.bsv-bs-popover-header
  // .fa-bullseye`) so the header is self-contained styling-wise.
  const popHeader = ui.divH([
    ui.iconFA('bullseye'),
    ui.span(['Binding Site']),
  ], 'bsv-bs-popover-header');
  popover.appendChild(popHeader);

  // Row: show binding site checkbox (the toggle).
  // `<label>` keeps the whole row clickable for toggling the checkbox without
  // needing `for=` / id plumbing — wraps both the text and the control cell.
  const showChk = document.createElement('input');
  showChk.type = 'checkbox';
  showChk.id = 'bsv-bs-show-chk';
  showChk.className = 'bsv-bs-show-chk';
  showChk.setAttribute('aria-label', 'Show side chains');
  const showRow = ui.element('label', 'bsv-bs-row') as HTMLLabelElement;
  showRow.appendChild(ui.divText('Show side chains', 'bsv-bs-row-label'));
  showRow.appendChild(ui.divH([showChk], 'bsv-bs-row-ctrl bsv-bs-row-ctrl--center'));

  // Row: radius slider. Raw `<input type=range>` kept deliberately — Datagrok's
  // `ui.input.slider` uses a Material blue accent + 28 px layout that clashes
  // with the Mol*-palette row we match to `.msp-control-row`.
  const radiusSlider = document.createElement('input');
  radiusSlider.type = 'range';
  radiusSlider.min = '3';
  radiusSlider.max = '10';
  radiusSlider.step = '0.5';
  radiusSlider.className = 'bsv-bs-radius-slider';
  radiusSlider.setAttribute('aria-label', 'Binding site radius in angstroms');
  const radiusVal = ui.divText('', 'bsv-bs-radius-val');
  const radiusRow = ui.divH([
    ui.divText('Radius', 'bsv-bs-row-label'),
    ui.divH([radiusSlider, radiusVal], 'bsv-bs-row-ctrl'),
  ], 'bsv-bs-row');

  const popBody = ui.divV([showRow, radiusRow]);
  popover.appendChild(popBody);

  // Popover is attached to the viewer root (same parent as wrapper) so it is
  // not affected by wrapper's pointer-events:none or flex layout. Attached
  // deferred in _bsvReposition, since wrapper.parentElement is only known
  // after the overlay is inserted into the DOM.

  /** Finds the visible bounding box for the popover — walks up from `parent`
   *  and intersects every overflow-hidden ancestor's rect. Mol*'s own
   *  layout wrappers (`.msp-layout-static`, `.msp-layout-region`) can be
   *  wider than the host Datagrok `<div>`, so taking only the first
   *  overflow:hidden ancestor would miss the tighter clip further up.
   *  @param {HTMLElement} parent element to start the walk from.
   *  @return {ClipBox} viewport-relative clip box (inclusive bounds). */
  type ClipBox = {left: number; right: number; top: number; bottom: number; width: number};
  const findClipBox = (parent: HTMLElement): ClipBox => {
    let minLeft = -Infinity;
    let minTop = -Infinity;
    let maxRight = Infinity;
    let maxBottom = Infinity;
    let el: HTMLElement | null = parent;
    while (el && el !== document.body) {
      const o = getComputedStyle(el);
      if (o.overflowX === 'hidden' || o.overflow === 'hidden') {
        const r = el.getBoundingClientRect();
        if (r.left > minLeft) minLeft = r.left;
        if (r.right < maxRight) maxRight = r.right;
      }
      if (o.overflowY === 'hidden' || o.overflow === 'hidden') {
        const r = el.getBoundingClientRect();
        if (r.top > minTop) minTop = r.top;
        if (r.bottom < maxBottom) maxBottom = r.bottom;
      }
      el = el.parentElement;
    }
    if (minLeft === -Infinity) {
      // No overflow-hidden ancestor — fall back to parent's rect.
      const r = parent.getBoundingClientRect();
      return {left: r.left, right: r.right, top: r.top, bottom: r.bottom, width: r.width};
    }
    return {left: minLeft, right: maxRight, top: minTop, bottom: maxBottom, width: maxRight - minLeft};
  };

  /** Position the popover in viewport coordinates (`position: fixed`). Body-
   *  attached to escape Mol*'s stacking context / overflow clip, but clamped
   *  to the `.msp-viewport` rectangle so it never floats outside the 3D
   *  panel's visible area. In a narrow Mol* viewport the popover shrinks;
   *  labels truncate with ellipsis, slider + value stay readable. */
  const positionPopover = () => {
    const parent = wrapper.parentElement;
    if (!parent) return;
    const iRect = iconBtn.getBoundingClientRect();
    // Viewport bounds: prefer `.msp-viewport`, fall back to the overlay's
    // own parent (which sits inside it).
    const viewport = (parent.closest('.msp-viewport') ?? parent) as HTMLElement;
    const vRect = viewport.getBoundingClientRect();

    const margin = 6;
    const naturalWidth = 290;
    const maxWidth = Math.max(160, vRect.width - margin * 2);
    const width = Math.min(naturalWidth, maxWidth);

    // Horizontal: anchor 36 px left of the icon, then clamp so the popover
    // stays inside [vRect.left + margin, vRect.right - margin].
    const preferredLeft = iRect.left - 36 - width;
    const minLeft = vRect.left + margin;
    const maxLeft = vRect.right - margin - width;
    const leftPx = Math.max(minLeft, Math.min(maxLeft, preferredLeft));

    // Vertical: pin to the TOP of the Mol* viewport — matches the native
    // `.msp-viewport-controls-panel` convention (Settings / Controls Info
    // panel). Using `strip.top` would place the popover at the icon's
    // vertical position, which is the *bottom* half of a tall narrow
    // viewport (the strip is bottom-anchored by Mol*'s layout).
    const topPx = vRect.top + margin;

    popover.style.width = `${width}px`;
    popover.style.left = `${leftPx}px`;
    popover.style.top = `${topPx}px`;
    popover.style.right = 'auto';
  };

  const updateSliderGradient = (r: number) => {
    // Only the percentage is a runtime value; the gradient + colours live in
    // CSS (`.bsv-bs-radius-slider`) via the `--bsv-bs-progress` custom property.
    const pct = ((r - 3) / 7 * 100).toFixed(1);
    radiusSlider.style.setProperty('--bsv-bs-progress', `${pct}%`);
  };

  const refreshUI = () => {
    const active = handlers.getShowBindingSite();
    const radius = handlers.getBindingSiteRadius();
    const ligAvail = handlers.isLigandAvailable();
    const applicable = handlers.hasPolymer();

    // No polymer loaded → hide the overlay entirely. Side-chain
    // highlighting is not meaningful without a protein.
    wrapper.classList.toggle('bsv-bs-overlay--hidden', !applicable);
    if (!applicable) {
      if (popoverOpen) closePopover();
      return;
    }

    // Swap Mol*'s native toggle-on/off classes so the icon inherits the
    // same color behavior as the native strip icons:
    //   off (inactive) → #9c835f   (muted — looks "transparent")
    //   on  (active)   → #332b1f   (fully saturated)
    //   :hover         → #ae5d04   (orange accent, automatic)
    iconBtn.classList.toggle('msp-btn-link-toggle-off', !active);
    iconBtn.classList.toggle('msp-btn-link-toggle-on', active);
    iconBtn.classList.toggle('bsv-bs-icon-btn--disabled', !ligAvail);
    iconBtn.disabled = !ligAvail;
    iconBtn.setAttribute('title', ligAvail ? 'Binding site' : 'No ligand detected');

    showChk.checked = active;
    showChk.disabled = !ligAvail;

    radiusSlider.value = String(radius);
    radiusVal.textContent = `${radius.toFixed(1)} \u00C5`;
    updateSliderGradient(radius);
  };

  let popoverOpen = false;

  const openPopover = () => {
    popoverOpen = true;
    // Attach to <body> so the popover escapes the Mol* viewport's stacking
    // context (and any overflow:hidden ancestor). Without this, in a narrow
    // docked viewer the popover gets clipped by the right-adjacent panel.
    // With position:fixed + body-attachment we can always render at full
    // 290 px width and overlap the strip / Components panel on top.
    if (popover.parentElement !== document.body) document.body.appendChild(popover);
    positionPopover();
    popover.classList.remove('bsv-bs-popover-hidden');
    refreshUI();
  };

  const closePopover = () => {
    popoverOpen = false;
    popover.classList.add('bsv-bs-popover-hidden');
  };

  // Icon click toggles the popover (no direct show/hide toggle on the icon).
  // Prevent focus-steal from Mol* canvas on mousedown.
  iconBtn.addEventListener('mousedown', (e: MouseEvent) => { e.preventDefault(); });
  iconBtn.addEventListener('click', (e: MouseEvent) => {
    e.stopPropagation();
    if (popoverOpen)
      closePopover();
    else
      openPopover();
    // Drop focus so Mol*'s generic `.msp-btn:focus` bg doesn't linger on the
    // icon (would leave the bullseye looking "pressed" until focus moves on).
    iconBtn.blur();
  });
  iconBtn.addEventListener('keydown', (e: KeyboardEvent) => {
    if (e.key === 'Enter' || e.key === ' ') {
      e.preventDefault();
      if (popoverOpen)
        closePopover();
      else
        openPopover();
    }
  });

  // Block mouse events from bubbling to the Mol* viewport so hover over our
  // overlay doesn't trigger residue highlighting underneath. Native Mol*
  // icons sit inside `.msp-viewport-controls` which Mol* special-cases; our
  // overlay is a sibling, so we must stop propagation ourselves.
  const swallowEvents = ['mousemove', 'pointermove', 'mouseover', 'mouseenter', 'wheel'];
  const swallow = (e: Event) => { e.stopPropagation(); };
  for (const evt of swallowEvents) {
    iconBtn.addEventListener(evt, swallow);
    popover.addEventListener(evt, swallow);
  }

  showChk.addEventListener('change', () => {
    if (!handlers.isLigandAvailable()) return;
    handlers.setShowBindingSite(showChk.checked);
    refreshUI();
  });

  let radiusTimer: ReturnType<typeof setTimeout> | null = null;
  radiusSlider.addEventListener('input', () => {
    const r = parseFloat(radiusSlider.value);
    radiusVal.textContent = `${r.toFixed(1)} \u00C5`;
    updateSliderGradient(r);
    if (radiusTimer !== null) clearTimeout(radiusTimer);
    radiusTimer = setTimeout(() => {
      handlers.setBindingSiteRadius(r);
      radiusTimer = null;
    }, 50);
  });

  const onKeydown = (e: KeyboardEvent) => {
    if (e.key === 'Escape' && popoverOpen) closePopover();
  };
  document.addEventListener('keydown', onKeydown, {capture: true});

  const onOutsideClick = (e: MouseEvent) => {
    if (!popoverOpen) return;
    const target = e.target as Node;
    if (wrapper.contains(target) || popover.contains(target)) return;
    closePopover();
  };
  document.addEventListener('click', onOutsideClick, {capture: false});

  // Re-measure whenever the viewer OR the Mol* strip resizes or moves. The
  // Mol* right sidebar can open/close without changing the viewer size but
  // it shifts the strip horizontally, so we observe the strip itself too.
  let resizeObserver: ResizeObserver | null = null;
  let mutationObserver: MutationObserver | null = null;
  let observedStrip: HTMLElement | null = null;
  let observedViewport: HTMLElement | null = null;
  let observedButtons: HTMLElement | null = null;
  const attachResizeObserver = () => {
    if (!wrapper.parentElement) return;
    if (typeof ResizeObserver !== 'undefined') {
      if (!resizeObserver) {
        resizeObserver = new ResizeObserver(() => positionOverlayRelativeToStrip(wrapper));
        resizeObserver.observe(wrapper.parentElement);
      }
      const strip = wrapper.parentElement.querySelector('.msp-viewport-controls') as HTMLElement | null;
      if (strip && strip !== observedStrip) {
        if (observedStrip) resizeObserver.unobserve(observedStrip);
        resizeObserver.observe(strip);
        observedStrip = strip;
      }
      const viewport = wrapper.parentElement.querySelector('.msp-viewport') as HTMLElement | null;
      if (viewport && viewport !== observedViewport) {
        if (observedViewport) resizeObserver.unobserve(observedViewport);
        resizeObserver.observe(viewport);
        observedViewport = viewport;
      }
    }
    // MutationObserver on the button stack catches late-rendered icons
    // (e.g. the selection-mode button). ResizeObserver doesn't always fire
    // before a late child's layout is measured.
    if (typeof MutationObserver !== 'undefined') {
      const buttons = wrapper.parentElement.querySelector('.msp-viewport-controls-buttons') as HTMLElement | null;
      if (buttons && buttons !== observedButtons) {
        if (mutationObserver) mutationObserver.disconnect();
        mutationObserver = new MutationObserver(() => positionOverlayRelativeToStrip(wrapper));
        mutationObserver.observe(buttons, {childList: true, subtree: true});
        observedButtons = buttons;
      }
    }
  };

  // Exposed so molstar-viewer.ts can trigger repositioning after the Mol*
  // plugin has finished rendering its chrome.
  const reposition = () => {
    positionOverlayRelativeToStrip(wrapper);
    attachResizeObserver();
  };
  wrapper._bsvReposition = () => {
    reposition();
    // Re-try a few times to catch async chrome mounts / panel animations.
    for (const ms of [150, 500, 1500]) setTimeout(reposition, ms);
  };

  // Safety net: a 500 ms polling interval. Wasteful but bulletproof against
  // any missed mutation/resize events. Stops when overlay is detached.
  const pollInterval = window.setInterval(() => {
    if (!wrapper.isConnected) {
      window.clearInterval(pollInterval);
      return;
    }
    positionOverlayRelativeToStrip(wrapper);
    attachResizeObserver();
    refreshUI();
  }, 500);

  wrapper._bsvCleanup = () => {
    document.removeEventListener('keydown', onKeydown, {capture: true});
    document.removeEventListener('click', onOutsideClick, {capture: false});
    if (radiusTimer !== null) clearTimeout(radiusTimer);
    if (resizeObserver) {
      resizeObserver.disconnect();
      resizeObserver = null;
    }
    if (mutationObserver) {
      mutationObserver.disconnect();
      mutationObserver = null;
    }
    if (popover.parentElement) popover.parentElement.removeChild(popover);
  };

  wrapper._bsvRefreshUI = refreshUI;

  refreshUI();

  return wrapper;
}
