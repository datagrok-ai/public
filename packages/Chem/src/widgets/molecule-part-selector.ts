/**
 * MoleculePartSelector — interactive 2D molecule view with atom selection.
 *
 * A reusable utility for picking parts of a molecule visually. Renders the
 * molecule via RDKit MinimalLib (SVG output), supports click-to-toggle and
 * box-drag selection, layered named highlights, and emits selection changes
 * through RxJS subjects compatible with existing Chem widgets.
 *
 * Designed to slot into the existing Datagrok chem pipeline:
 *   - emits / consumes the same {@link ISubstruct} shape used by the grid
 *     cell renderer
 *   - composed selection + highlight layers can be exposed as an
 *     {@link ISubstructProvider}, picked up automatically by the chem
 *     `RDKitCellRenderer` when registered on a column's temp store
 *   - uses the shared `getRdKitModule()` so it shares RDKit state with the
 *     rest of the package
 *
 * Typical use:
 * ```ts
 * import {createMoleculePartSelector} from './widgets/molecule-part-selector';
 *
 * const selector = createMoleculePartSelector({
 *   molecule: 'CC(=O)Oc1ccccc1C(=O)O',
 *   width: 320,
 *   height: 240,
 *   mode: 'both',
 * });
 * container.appendChild(selector.root);
 *
 * selector.onSelectionChanged.subscribe((atoms) => {
 *   console.log('selected', atoms);
 * });
 *
 * selector.highlight({
 *   id: 'pka',
 *   atoms: [16, 23],
 *   color: 'rgba(255, 153, 51, 1)',
 *   notes: {16: '2.64', 23: '9.90'},
 * });
 * ```
 */

import {BehaviorSubject, Subject} from 'rxjs';

import {ISubstruct, ISubstructProvider} from '@datagrok-libraries/chem-meta/src/types';
import {RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';

import {getRdKitModule} from '../utils/chem-common-rdkit';
import {hexToPercentRgb} from '../utils/chem-common';

/**
 * Selection input mode.
 *
 * - `'click'`  — click an atom to toggle it
 * - `'box'`    — drag a rectangle to select all atoms inside it
 * - `'lasso'`  — trace a free-form polygon to select atoms inside it
 * - `'both'`   — click + box-drag (lasso must be enabled explicitly)
 * - `'all'`    — click + box-drag + lasso (modifier key chooses between box and lasso)
 * - `'none'`   — render only, no input
 *
 * In `'all'` mode, holding `Ctrl`/`Cmd` while dragging starts a lasso instead
 * of a rectangle. Modifier semantics for both box and lasso are identical:
 * no modifier = replace, `Shift` = add, `Alt` = subtract.
 */
export type SelectionMode = 'click' | 'box' | 'lasso' | 'both' | 'all' | 'none';

/** A named highlight layer composed into the rendered molecule. */
export interface HighlightLayer {
  /** Stable identifier — re-adding a layer with the same id replaces it. */
  id: string;
  /** Atoms to highlight. */
  atoms: number[];
  /** Optional bonds to highlight (in addition to bonds inferred between atoms). */
  bonds?: number[];
  /** Color as hex string ('#ff9933'), 'rgba(...)', or as a [r,g,b,a] array in 0..1. */
  color: number[] | string;
  /** Per-atom annotation strings rendered next to atoms (e.g. pKa, R1, score). */
  notes?: {[atomIdx: number]: string};
  /** Per-bond annotation strings. */
  bondNotes?: {[bondIdx: number]: string};
  /** Higher priority wins when atoms appear in multiple layers (default 0). */
  priority?: number;
}

/** Constructor options for {@link MoleculePartSelector}. */
export interface MoleculePartSelectorOptions {
  /** SMILES or molblock to render. */
  molecule: string;
  /** Initial host width in pixels. Default 320. */
  width?: number;
  /** Initial host height in pixels. Default 240. */
  height?: number;
  /** Selection mode. Default `'both'`. */
  mode?: SelectionMode;
  /** Atoms to pre-select. */
  initialSelection?: number[];
  /** Layers to add at construction time. */
  initialLayers?: HighlightLayer[];
  /** Color of the user selection (default magenta). */
  selectionColor?: number[] | string;
  /** Whether to overlay {@link HighlightLayer.notes} as text labels. Default true. */
  showNotes?: boolean;
}

const DEFAULT_SELECTION_COLOR: number[] = [0.85, 0.4, 0.75, 1.0];

/** Normalize a color from hex / rgba string / number array to [r,g,b,a] in 0..1. */
function normalizeColor(c: number[] | string | undefined, fallback: number[]): number[] {
  if (!c) return fallback;
  if (Array.isArray(c)) return c;
  if (typeof c === 'string') {
    if (c.startsWith('#')) return hexToPercentRgb(c) ?? fallback;
    const m = /rgba?\s*\(\s*([\d.]+)\s*,\s*([\d.]+)\s*,\s*([\d.]+)(?:\s*,\s*([\d.]+))?\s*\)/i.exec(c);
    if (m) {
      return [
        parseFloat(m[1]) / 255,
        parseFloat(m[2]) / 255,
        parseFloat(m[3]) / 255,
        m[4] !== undefined ? parseFloat(m[4]) : 1.0,
      ];
    }
  }
  return fallback;
}

/**
 * Interactive molecule view with atom selection. See module-level doc above.
 */
export class MoleculePartSelector {
  /** The DOM element to insert into the page. */
  readonly root: HTMLElement;

  /** Fires whenever the selection changes (init value: empty array). */
  readonly onSelectionChanged: BehaviorSubject<number[]>;
  /** Fires once per atom click, regardless of whether the click added or removed it. */
  readonly onAtomClicked: Subject<number>;
  /** Fires when the rendered molecule is replaced via {@link setMolecule}. */
  readonly onMoleculeChanged: Subject<string>;

  // ---- private state -------------------------------------------------------

  private _molecule: string = '';
  private _mol: RDMol | null = null;
  private readonly _layers: Map<string, HighlightLayer> = new Map();
  private _selection: Set<number> = new Set();
  private _mode: SelectionMode;
  private readonly _width: number;
  private readonly _height: number;
  private readonly _selectionColor: number[];
  private readonly _showNotes: boolean;
  private _disposed = false;
  private _renderPending = false;
  private _resizeObserver: ResizeObserver | null = null;

  // box-drag and lasso state
  private _dragStart: {x: number; y: number} | null = null;
  private _dragRect: HTMLDivElement | null = null;
  private _dragModifiers: {shift: boolean; alt: boolean} | null = null;
  private _suppressNextClick = false;

  // lasso-specific
  private _lassoActive = false;
  private _lassoPoints: Array<{x: number; y: number}> = [];
  private _lassoSvg: SVGSVGElement | null = null;
  private _lassoPath: SVGPathElement | null = null;
  private static readonly _LASSO_MIN_STEP = 4; // px between sampled points

  // bound listener references for clean removal in dispose()
  private readonly _onMouseDown: (e: MouseEvent) => void;
  private readonly _onMouseMove: (e: MouseEvent) => void;
  private readonly _onMouseUp: (e: MouseEvent) => void;
  private readonly _onClick: (e: MouseEvent) => void;

  constructor(opts: MoleculePartSelectorOptions) {
    this._width = opts.width ?? 320;
    this._height = opts.height ?? 240;
    this._mode = opts.mode ?? 'both';
    this._selectionColor = normalizeColor(opts.selectionColor, DEFAULT_SELECTION_COLOR);
    this._showNotes = opts.showNotes ?? true;

    this.root = document.createElement('div');
    this.root.className = 'chem-mol-part-selector';
    Object.assign(this.root.style, {
      position: 'relative',
      width: this._width + 'px',
      height: this._height + 'px',
      userSelect: 'none',
      WebkitUserSelect: 'none',
    } as Partial<CSSStyleDeclaration>);

    this.onSelectionChanged = new BehaviorSubject<number[]>([]);
    this.onAtomClicked = new Subject<number>();
    this.onMoleculeChanged = new Subject<string>();

    this._onMouseDown = (e) => this._handleMouseDown(e);
    this._onMouseMove = (e) => this._handleMouseMove(e);
    this._onMouseUp = (e) => this._handleMouseUp(e);
    this._onClick = (e) => this._handleClick(e);

    this.root.addEventListener('mousedown', this._onMouseDown);
    document.addEventListener('mousemove', this._onMouseMove);
    document.addEventListener('mouseup', this._onMouseUp);
    this.root.addEventListener('click', this._onClick);

    if (opts.initialLayers) {
      for (const layer of opts.initialLayers)
        this._layers.set(layer.id, layer);
    }
    if (opts.initialSelection) {
      for (const a of opts.initialSelection)
        this._selection.add(a);
    }

    if (typeof ResizeObserver !== 'undefined') {
      this._resizeObserver = new ResizeObserver(() => this._scheduleRender());
      this._resizeObserver.observe(this.root);
    }

    this.setMolecule(opts.molecule);
  }

  // ---- public API: molecule ------------------------------------------------

  /** Replaces the rendered molecule. Accepts SMILES or molblock. */
  setMolecule(molecule: string): void {
    if (this._disposed) return;
    if (this._mol) {
      try { this._mol.delete(); } catch { /* ignore */ }
      this._mol = null;
    }
    this._molecule = molecule;
    try {
      const mod = getRdKitModule();
      const m = mod.get_mol(molecule);
      if (m && m.is_valid()) {
        this._mol = m;
      } else {
        m && m.delete();
        this._mol = null;
      }
    } catch {
      this._mol = null;
    }
    this.onMoleculeChanged.next(molecule);
    this._scheduleRender();
  }

  /** Returns the currently rendered molecule string (as last passed in). */
  getMolecule(): string {
    return this._molecule;
  }

  // ---- public API: selection ----------------------------------------------

  /** Sets the selection mode. */
  setMode(mode: SelectionMode): void {
    this._mode = mode;
  }

  /** Returns the current selection mode. */
  getMode(): SelectionMode {
    return this._mode;
  }

  /** Replaces the entire selection. */
  setSelection(atoms: number[]): void {
    this._selection = new Set(atoms);
    this._notifySelection();
    this._scheduleRender();
  }

  /** Adds atoms to the selection (no-op for atoms already selected). */
  addToSelection(atoms: number[]): void {
    let changed = false;
    for (const a of atoms) {
      if (!this._selection.has(a)) { this._selection.add(a); changed = true; }
    }
    if (changed) {
      this._notifySelection();
      this._scheduleRender();
    }
  }

  /** Removes atoms from the selection. */
  removeFromSelection(atoms: number[]): void {
    let changed = false;
    for (const a of atoms) {
      if (this._selection.delete(a)) changed = true;
    }
    if (changed) {
      this._notifySelection();
      this._scheduleRender();
    }
  }

  /** Returns the current selection as a sorted array of atom indices. */
  getSelection(): number[] {
    return [...this._selection].sort((a, b) => a - b);
  }

  /** Clears the selection. Does not affect highlight layers. */
  clear(): void {
    if (this._selection.size === 0) return;
    this._selection.clear();
    this._notifySelection();
    this._scheduleRender();
  }

  // ---- public API: highlight layers ---------------------------------------

  /** Adds or replaces a named highlight layer. */
  highlight(layer: HighlightLayer): void {
    this._layers.set(layer.id, layer);
    this._scheduleRender();
  }

  /** Removes a highlight layer by id. */
  removeHighlight(id: string): void {
    if (this._layers.delete(id))
      this._scheduleRender();
  }

  /** Removes all highlight layers. Does not affect the selection. */
  clearHighlights(): void {
    if (this._layers.size === 0) return;
    this._layers.clear();
    this._scheduleRender();
  }

  /** Returns a snapshot of all current highlight layers. */
  getHighlights(): HighlightLayer[] {
    return [...this._layers.values()];
  }

  // ---- public API: export -------------------------------------------------

  /**
   * Returns the current selection as a stand-alone {@link ISubstruct}, ready
   * to feed back into the cell renderer / providers system.
   */
  getSelectionAsSubstruct(): ISubstruct {
    const atoms = this.getSelection();
    const highlightAtomColors: {[k: number]: number[]} = {};
    for (const a of atoms) highlightAtomColors[a] = this._selectionColor;
    return {atoms, highlightAtomColors};
  }

  /**
   * Returns selection + all visible highlight layers composed into a single
   * {@link ISubstruct}. This is exactly what the SVG renderer below feeds to
   * RDKit, so consumers (e.g. an `ISubstructProvider` registered on a column)
   * can produce identical highlights in the grid.
   */
  getCombinedSubstruct(): ISubstruct {
    return this._composeSubstruct();
  }

  /**
   * Wraps this selector as an {@link ISubstructProvider} that returns the
   * combined substruct for a single row. Use with `addSubstructProvider` from
   * `@datagrok-libraries/chem-meta/src/types` to mirror the picker's
   * highlights into a Datagrok grid column.
   *
   * @param targetRowIdx if provided, only that row receives the substruct.
   */
  asSubstructProvider(targetRowIdx: number | null = null): ISubstructProvider {
    return {
      getSubstruct: (rowIdx: number | null): ISubstruct | undefined => {
        if (targetRowIdx !== null && rowIdx !== targetRowIdx) return undefined;
        return this._composeSubstruct();
      },
    };
  }

  // ---- public API: lifecycle ----------------------------------------------

  /**
   * Frees the underlying RDMol, removes DOM listeners, and completes all
   * observable subjects. The instance must not be used after dispose().
   */
  dispose(): void {
    if (this._disposed) return;
    this._disposed = true;
    this.root.removeEventListener('mousedown', this._onMouseDown);
    document.removeEventListener('mousemove', this._onMouseMove);
    document.removeEventListener('mouseup', this._onMouseUp);
    this.root.removeEventListener('click', this._onClick);
    if (this._resizeObserver) {
      this._resizeObserver.disconnect();
      this._resizeObserver = null;
    }
    if (this._mol) {
      try { this._mol.delete(); } catch { /* ignore */ }
      this._mol = null;
    }
    this.onSelectionChanged.complete();
    this.onAtomClicked.complete();
    this.onMoleculeChanged.complete();
    this.root.innerHTML = '';
  }

  // ---- private: substructure composition ----------------------------------

  private _notifySelection(): void {
    this.onSelectionChanged.next(this.getSelection());
  }

  private _composeSubstruct(): ISubstruct {
    const atoms: number[] = [];
    const bonds: number[] = [];
    const highlightAtomColors: {[k: number]: number[]} = {};
    const highlightBondColors: {[k: number]: number[]} = {};
    const atomNotes: {[k: number]: string} = {};
    const bondNotes: {[k: number]: string} = {};

    // sort layers ascending by priority so higher priority overrides on overlap
    const sorted = [...this._layers.values()].sort(
      (a, b) => (a.priority ?? 0) - (b.priority ?? 0));
    for (const layer of sorted) {
      const color = normalizeColor(layer.color, this._selectionColor);
      for (const a of layer.atoms) {
        atoms.push(a);
        highlightAtomColors[a] = color;
      }
      if (layer.bonds) {
        for (const b of layer.bonds) {
          bonds.push(b);
          highlightBondColors[b] = color;
        }
      }
      if (layer.notes) Object.assign(atomNotes, layer.notes);
      if (layer.bondNotes) Object.assign(bondNotes, layer.bondNotes);
    }

    // user selection always wins over layers
    for (const a of this._selection) {
      atoms.push(a);
      highlightAtomColors[a] = this._selectionColor;
    }

    return {atoms, bonds, highlightAtomColors, highlightBondColors, atomNotes, bondNotes};
  }

  // ---- private: rendering -------------------------------------------------

  private _scheduleRender(): void {
    if (this._disposed || this._renderPending) return;
    this._renderPending = true;
    requestAnimationFrame(() => {
      this._renderPending = false;
      if (!this._disposed) this._render();
    });
  }

  private _render(): void {
    if (!this._mol || !this._mol.is_valid()) {
      this.root.innerHTML = '<div style="color:#a00;padding:8px;font-size:11px">' +
        'Failed to render molecule</div>';
      return;
    }
    const sub = this._composeSubstruct();
    const w = this.root.clientWidth || this._width;
    const h = this.root.clientHeight || this._height;

    const details: {[k: string]: any} = {
      width: w,
      height: h,
      atoms: sub.atoms,
      bonds: sub.bonds,
      highlightAtomColors: sub.highlightAtomColors,
      highlightBondColors: sub.highlightBondColors,
      bondLineWidth: 1.0,
      multipleBondOffset: 0.18,
      baseFontSize: 0.6,
    };

    let svgString = '';
    try {
      svgString = this._mol.get_svg_with_highlights(JSON.stringify(details));
    } catch {
      this.root.innerHTML = '<div style="color:#a00;padding:8px;font-size:11px">' +
        'SVG render failed</div>';
      return;
    }

    this.root.innerHTML = svgString;
    const svgEl = this.root.querySelector('svg') as SVGSVGElement | null;
    if (!svgEl) return;
    svgEl.style.display = 'block';
    svgEl.style.width = '100%';
    svgEl.style.height = '100%';

    if (this._showNotes && sub.atomNotes && Object.keys(sub.atomNotes).length > 0)
      this._overlayNotes(svgEl, sub);
  }

  private _overlayNotes(svgEl: SVGSVGElement, sub: ISubstruct): void {
    const positions = this._getHostAtomPositions(svgEl);
    for (const [idxStr, value] of Object.entries(sub.atomNotes ?? {})) {
      const idx = parseInt(idxStr, 10);
      const p = positions.get(idx);
      if (!p) continue;
      const note = document.createElement('div');
      note.className = 'chem-mps-note';
      note.textContent = String(value);
      Object.assign(note.style, {
        position: 'absolute',
        left: p.x + 'px',
        top: p.y + 'px',
        transform: 'translate(-50%, -110%)',
        font: '600 10px monospace',
        background: 'rgba(255,255,255,0.88)',
        padding: '1px 3px',
        borderRadius: '3px',
        pointerEvents: 'none',
        whiteSpace: 'nowrap',
        color: '#1d2330',
      } as Partial<CSSStyleDeclaration>);
      this.root.appendChild(note);
    }
  }

  // ---- private: SVG -> atom-position mapping ------------------------------

  private _extractSvgAtomPositions(svgEl: SVGSVGElement): Map<number, {x: number; y: number}> {
    const positions = new Map<number, {x: number; y: number}>();

    // 1) highlighted atoms — RDKit emits <ellipse class="atom-N" cx="..." cy="..."/>
    //    one per highlighted atom. This is the most accurate signal because
    //    the ellipse is anchored exactly at the atom centre (vs <text> which
    //    gives a bbox center, vs bond endpoints which are slightly shrunk).
    const ellipses = svgEl.querySelectorAll('ellipse[class*="atom-"]');
    for (let i = 0; i < ellipses.length; i++) {
      const el = ellipses[i];
      const cls = el.getAttribute('class') || '';
      // skip if the class also contains a bond- token (defensive — shouldn't happen)
      if (/(?:^|\s)bond-\d+/.test(cls)) continue;
      const m = /(?:^|\s)atom-(\d+)/.exec(cls);
      if (!m) continue;
      const idx = parseInt(m[1], 10);
      const cx = parseFloat(el.getAttribute('cx') || 'NaN');
      const cy = parseFloat(el.getAttribute('cy') || 'NaN');
      if (!Number.isNaN(cx) && !Number.isNaN(cy))
        positions.set(idx, {x: cx, y: cy});
    }

    // 2) heteroatoms with text labels — RDKit emits <text class="atom-N">
    const texts = svgEl.querySelectorAll('text[class*="atom-"]');
    for (let i = 0; i < texts.length; i++) {
      const t = texts[i];
      const cls = t.getAttribute('class') || '';
      const m = /(?:^|\s)atom-(\d+)/.exec(cls);
      if (!m) continue;
      const idx = parseInt(m[1], 10);
      if (positions.has(idx)) continue; // ellipse wins (more accurate)
      try {
        const bb = (t as unknown as SVGGraphicsElement).getBBox();
        positions.set(idx, {x: bb.x + bb.width / 2, y: bb.y + bb.height / 2});
      } catch { /* ignore */ }
    }

    // 3) carbons / unlabeled atoms — average bond endpoints from <path class="bond-K atom-A atom-B">.
    //    Used as the final fallback for atoms with neither an ellipse highlight
    //    nor a text label.
    const bondEnds = new Map<number, Array<{x: number; y: number}>>();
    const bondPaths = svgEl.querySelectorAll('path[class*="bond-"]');
    for (let i = 0; i < bondPaths.length; i++) {
      const p = bondPaths[i];
      const cls = p.getAttribute('class') || '';
      const ids: number[] = [];
      const matches = cls.match(/atom-(\d+)/g) || [];
      for (const mm of matches) {
        const n = parseInt(mm.slice(5), 10);
        if (!Number.isNaN(n)) ids.push(n);
      }
      if (ids.length !== 2) continue;
      const d = p.getAttribute('d') || '';
      const coords: Array<{x: number; y: number}> = [];
      const re = /[ML]\s*([\-\d.]+)[\s,]+([\-\d.]+)/g;
      let mm: RegExpExecArray | null;
      while ((mm = re.exec(d)) !== null)
        coords.push({x: parseFloat(mm[1]), y: parseFloat(mm[2])});
      if (coords.length < 2) continue;
      if (!bondEnds.has(ids[0])) bondEnds.set(ids[0], []);
      if (!bondEnds.has(ids[1])) bondEnds.set(ids[1], []);
      bondEnds.get(ids[0])!.push(coords[0]);
      bondEnds.get(ids[1])!.push(coords[coords.length - 1]);
    }
    for (const [idx, pts] of bondEnds.entries()) {
      if (positions.has(idx)) continue;
      const cx = pts.reduce((s, p) => s + p.x, 0) / pts.length;
      const cy = pts.reduce((s, p) => s + p.y, 0) / pts.length;
      positions.set(idx, {x: cx, y: cy});
    }

    return positions;
  }

  private _getHostAtomPositions(svgEl: SVGSVGElement): Map<number, {x: number; y: number}> {
    const svgPositions = this._extractSvgAtomPositions(svgEl);
    const vb = svgEl.viewBox && svgEl.viewBox.baseVal;
    const svgRect = svgEl.getBoundingClientRect();
    const hostRect = this.root.getBoundingClientRect();
    const sx = vb && vb.width ? svgRect.width / vb.width : 1;
    const sy = vb && vb.height ? svgRect.height / vb.height : 1;
    const dx = (svgRect.left - hostRect.left) - (vb ? vb.x * sx : 0);
    const dy = (svgRect.top - hostRect.top) - (vb ? vb.y * sy : 0);
    const out = new Map<number, {x: number; y: number}>();
    for (const [idx, p] of svgPositions.entries())
      out.set(idx, {x: dx + p.x * sx, y: dy + p.y * sy});
    return out;
  }

  private _findAtomFromTarget(target: EventTarget | null): number | null {
    let el: any = target;
    while (el && el !== this.root && el.tagName !== 'svg') {
      const cls = el.getAttribute && el.getAttribute('class');
      const m = cls && /(?:^|\s)atom-(\d+)/.exec(cls);
      if (m) return parseInt(m[1], 10);
      el = el.parentNode;
    }
    return null;
  }

  // ---- private: input handlers --------------------------------------------

  private _handleClick(e: MouseEvent): void {
    if (this._suppressNextClick) {
      this._suppressNextClick = false;
      return;
    }
    if (this._mode !== 'click' && this._mode !== 'both') return;
    const idx = this._findAtomFromTarget(e.target);
    if (idx === null) return;
    if (this._selection.has(idx)) this._selection.delete(idx);
    else this._selection.add(idx);
    this.onAtomClicked.next(idx);
    this._notifySelection();
    this._scheduleRender();
  }

  /** Decide whether the current mode + modifiers should start a drag at all,
   *  and if so, whether it's a box-drag or a lasso. */
  private _dragKindForMouseDown(e: MouseEvent): 'box' | 'lasso' | null {
    const ctrlOrCmd = e.ctrlKey || e.metaKey;
    switch (this._mode) {
      case 'box':   return 'box';
      case 'lasso': return 'lasso';
      case 'both':  return 'box';
      case 'all':   return ctrlOrCmd ? 'lasso' : 'box';
      default:      return null; // 'click' / 'none'
    }
  }

  private _handleMouseDown(e: MouseEvent): void {
    if (e.button !== 0) return;
    const kind = this._dragKindForMouseDown(e);
    if (kind === null) return;
    // let single-atom clicks fall through to the click handler
    if (this._findAtomFromTarget(e.target) !== null) return;

    const hostRect = this.root.getBoundingClientRect();
    const start = {x: e.clientX - hostRect.left, y: e.clientY - hostRect.top};
    this._dragStart = start;
    this._dragModifiers = {shift: e.shiftKey, alt: e.altKey};

    if (kind === 'box') {
      const rect = document.createElement('div');
      rect.className = 'chem-mps-drag-rect';
      Object.assign(rect.style, {
        position: 'absolute',
        left: start.x + 'px',
        top: start.y + 'px',
        width: '0px',
        height: '0px',
        border: '1px dashed #4a8',
        background: 'rgba(102, 204, 102, 0.12)',
        pointerEvents: 'none',
        boxSizing: 'border-box',
      } as Partial<CSSStyleDeclaration>);
      this.root.appendChild(rect);
      this._dragRect = rect;
    } else {
      this._lassoActive = true;
      this._lassoPoints = [start];
      this._createLassoOverlay();
    }
    e.preventDefault();
  }

  private _handleMouseMove(e: MouseEvent): void {
    if (!this._dragStart) return;
    const hostRect = this.root.getBoundingClientRect();
    const cx = e.clientX - hostRect.left;
    const cy = e.clientY - hostRect.top;

    if (this._lassoActive) {
      const last = this._lassoPoints[this._lassoPoints.length - 1];
      // throttle: only sample when the cursor has moved a few px
      if (Math.hypot(cx - last.x, cy - last.y) >= MoleculePartSelector._LASSO_MIN_STEP) {
        this._lassoPoints.push({x: cx, y: cy});
        this._updateLassoOverlay();
      }
      return;
    }

    if (!this._dragRect) return;
    const x = Math.min(this._dragStart.x, cx);
    const y = Math.min(this._dragStart.y, cy);
    this._dragRect.style.left = x + 'px';
    this._dragRect.style.top = y + 'px';
    this._dragRect.style.width = Math.abs(cx - this._dragStart.x) + 'px';
    this._dragRect.style.height = Math.abs(cy - this._dragStart.y) + 'px';
  }

  private _handleMouseUp(e: MouseEvent): void {
    if (!this._dragStart) return;
    const hostRect = this.root.getBoundingClientRect();
    const cx = e.clientX - hostRect.left;
    const cy = e.clientY - hostRect.top;

    let moved = false;

    if (this._lassoActive) {
      this._lassoPoints.push({x: cx, y: cy});
      // need at least a small triangle to enclose anything
      moved = this._lassoPathLength() > 12 && this._lassoPoints.length >= 3;
      if (moved) this._applyLassoSelection();
    } else if (this._dragRect) {
      const x0 = Math.min(this._dragStart.x, cx);
      const x1 = Math.max(this._dragStart.x, cx);
      const y0 = Math.min(this._dragStart.y, cy);
      const y1 = Math.max(this._dragStart.y, cy);
      moved = Math.hypot(x1 - x0, y1 - y0) > 4;
      if (moved) this._applyBoxSelection(x0, y0, x1, y1);
    }

    if (moved) {
      // suppress the synthetic click that follows the drag mouseup
      this._suppressNextClick = true;
      setTimeout(() => { this._suppressNextClick = false; }, 0);
    }

    this._dragStart = null;
    this._dragModifiers = null;
    if (this._dragRect && this._dragRect.parentNode)
      this._dragRect.parentNode.removeChild(this._dragRect);
    this._dragRect = null;
    this._lassoActive = false;
    this._lassoPoints = [];
    this._destroyLassoOverlay();

    if (moved) {
      this._notifySelection();
      this._scheduleRender();
    }
  }

  // ---- private: box and lasso application --------------------------------

  private _applyBoxSelection(x0: number, y0: number, x1: number, y1: number): void {
    const svgEl = this.root.querySelector('svg') as SVGSVGElement | null;
    if (!svgEl) return;
    const positions = this._getHostAtomPositions(svgEl);
    if (!this._dragModifiers!.shift && !this._dragModifiers!.alt)
      this._selection.clear();
    for (const [idx, p] of positions.entries()) {
      if (p.x >= x0 && p.x <= x1 && p.y >= y0 && p.y <= y1) {
        if (this._dragModifiers!.alt) this._selection.delete(idx);
        else this._selection.add(idx);
      }
    }
  }

  private _applyLassoSelection(): void {
    const svgEl = this.root.querySelector('svg') as SVGSVGElement | null;
    if (!svgEl) return;
    const positions = this._getHostAtomPositions(svgEl);
    const poly = this._lassoPoints;
    if (!this._dragModifiers!.shift && !this._dragModifiers!.alt)
      this._selection.clear();
    for (const [idx, p] of positions.entries()) {
      if (MoleculePartSelector._pointInPolygon(p.x, p.y, poly)) {
        if (this._dragModifiers!.alt) this._selection.delete(idx);
        else this._selection.add(idx);
      }
    }
  }

  /** Standard ray-casting point-in-polygon test. Works for any simple polygon
   *  including non-convex shapes. O(n) per point against an n-vertex polygon. */
  private static _pointInPolygon(
    x: number, y: number, poly: Array<{x: number; y: number}>): boolean {
    let inside = false;
    for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
      const xi = poly[i].x; const yi = poly[i].y;
      const xj = poly[j].x; const yj = poly[j].y;
      if (((yi > y) !== (yj > y)) &&
          (x < (xj - xi) * (y - yi) / (yj - yi) + xi))
        inside = !inside;
    }
    return inside;
  }

  private _lassoPathLength(): number {
    let len = 0;
    for (let i = 1; i < this._lassoPoints.length; i++) {
      const a = this._lassoPoints[i - 1];
      const b = this._lassoPoints[i];
      len += Math.hypot(b.x - a.x, b.y - a.y);
    }
    return len;
  }

  // ---- private: lasso overlay ---------------------------------------------

  private _createLassoOverlay(): void {
    // Use a transparent SVG positioned over the host to draw the lasso path.
    // Plain HTML can't render arbitrary polylines easily; SVG is the right tool.
    const SVG_NS = 'http://www.w3.org/2000/svg';
    const overlay = document.createElementNS(SVG_NS, 'svg');
    overlay.setAttribute('class', 'chem-mps-lasso-overlay');
    Object.assign(overlay.style, {
      position: 'absolute',
      left: '0',
      top: '0',
      width: '100%',
      height: '100%',
      pointerEvents: 'none',
    } as Partial<CSSStyleDeclaration>);
    const path = document.createElementNS(SVG_NS, 'path');
    path.setAttribute('fill', 'rgba(102, 204, 102, 0.12)');
    path.setAttribute('stroke', '#4a8');
    path.setAttribute('stroke-width', '1.2');
    path.setAttribute('stroke-dasharray', '4 3');
    path.setAttribute('d', '');
    overlay.appendChild(path);
    this.root.appendChild(overlay);
    this._lassoSvg = overlay as unknown as SVGSVGElement;
    this._lassoPath = path;
    this._updateLassoOverlay();
  }

  private _updateLassoOverlay(): void {
    if (!this._lassoPath || this._lassoPoints.length === 0) return;
    let d = 'M ' + this._lassoPoints[0].x + ' ' + this._lassoPoints[0].y;
    for (let i = 1; i < this._lassoPoints.length; i++)
      d += ' L ' + this._lassoPoints[i].x + ' ' + this._lassoPoints[i].y;
    if (this._lassoPoints.length > 2) d += ' Z';
    this._lassoPath.setAttribute('d', d);
  }

  private _destroyLassoOverlay(): void {
    if (this._lassoSvg && this._lassoSvg.parentNode)
      this._lassoSvg.parentNode.removeChild(this._lassoSvg);
    this._lassoSvg = null;
    this._lassoPath = null;
  }
}

/**
 * Convenience factory — equivalent to `new MoleculePartSelector(opts)`.
 *
 * Provided for ergonomics; both forms are supported.
 */
export function createMoleculePartSelector(
  opts: MoleculePartSelectorOptions): MoleculePartSelector {
  return new MoleculePartSelector(opts);
}
