import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Subscription} from 'rxjs';

import {drawMoleculeToCanvas} from '../../utils/chem-common-rdkit';
import {SarMatrix} from './sar-matrix-types';

/** Layout, colour and grid types shared by the viewer and its panels. */

/** Transparent (alpha 0) so a drawn core blends with the card/pane instead of showing a white box. */
export const CORE_BG_ARGB = 0x00000000;
export const HEADER_ARGB = 0xFFF7F7F9;
export const WHITE_ARGB = 0xFFFFFFFF;
export const ACTIVITY_SCHEME = [DG.Color.red, DG.Color.green];

export const CELL_W = 104;
export const CELL_H = 76;
export const CELL_W_MAX = 210;
export const CHIP_H = 13;
export const CHIP_PAD = 3;
export const CHIP_MARGIN = 3;
export const HEADER_W = 78;
export const HEADER_H = 46;
export const COL_HEADER_H = HEADER_H + 36;
export const CORE_W = 132;
/** Room for the grid's own horizontal scrollbar when a grid is sized to its rows, not stretched. */
export const GRID_SCROLLBAR_H = 18;
export const BENEFIT_MOL_W = 62;
export const BENEFIT_MOL_H = 34;
export const CARD_CORE_W = 78;
export const CARD_CORE_H = 44;
/** A whole assembled analog needs more room than the core or substituent it was built from. */
export const ANALOG_W = 220;
/** Must track the `.chem-sar-nav` width in the stylesheet, or cells are fitted against the wrong pane. */
export const NAV_W = 320;
export const NAV_COLLAPSED_W = 28;
export const TABLE_CHROME = 60;

/** Cell-tint alpha (0-255): solid for observed, fainter for virtual, faintest for thin predictions. */
export const REAL_ALPHA = 102;
export const VIRTUAL_ALPHA = 46;
export const VIRTUAL_ALPHA_MIN = 16;
/** Support at/above which a virtual cell is at full alpha. */
export const FULL_SUPPORT = 3;

export const TAB_MATRIX = 'SAR Matrix';
export const TAB_TRANSFER = 'SAR Transfer';
export const TAB_MAKELIST = 'Make list';

/** A cell addressed by the matrix it belongs to and its position in it — how the panels hand single
 *  cells to one another, since a cell object alone cannot say where it came from. */
export interface MatrixCellRef {
  matrix: SarMatrix;
  ri: number;
  ci: number;
}

/** One displayed grid row. Each row carries its OWN matrix, because a transfer's two sides can come
 *  from different matrices. */
export interface PaneRow {
  matrix: SarMatrix;
  rowIndex: number;
  colIdxs: number[];
  label: string;
  /** Highlight the row's most potent observed cell; set for a transfer's two sides. */
  markBest?: boolean;
}

export interface PaneColumn {
  substSmiles: string;
  position: string;
  /** Metric caption under the depiction; empty when the Label control is None. */
  caption: string;
}

/** The rendered pane grid plus per-render row/column state, computed once at build time since every
 *  paint reads it. */
export interface MatrixGridState {
  grid: DG.Grid;
  /** Scaffold frame backing the grid; held so the grid keeps a live reference for its lifetime. */
  df: DG.DataFrame;
  rows: PaneRow[];
  columns: PaneColumn[];
  /** Grid column NAME (never the grid column index, which pinning shifts) -> index into `columns`. */
  colKeyToIdx: Map<string, number>;
  /** Indices into `columns` that begin an R-position group — only those draw the position label. */
  firstOfGroup: Set<number>;
  /** What the Core header draws: the shared core where the rows have one, else the matrix's own key,
   *  which is the usual case since a matrix's rows are different cores by construction. */
  headerDepiction: string | null;
  /** One alignment template for the whole pane so an attachment point sits in the same place
   *  everywhere. Null when rows span several matrices and share no core. */
  paneTemplate: string | null;
  /** Per displayed row, the molblock its cells align to; filled on first paint and kept for the
   *  life of the grid. */
  rowTemplates: (string | null)[];
}

/** A pane's grid plus the subscriptions that must be torn down with it. */
export interface PaneGridSlot {
  state: MatrixGridState | null;
  subs: {unsubscribe(): void}[];
}

export const GRID_FONT = 'Roboto, "Segoe UI", sans-serif';

/** Resolve a Datagrok CSS palette variable to a color string, memoized so per-cell paints don't
 *  re-run `getComputedStyle`. */
const cssColorCache = new Map<string, string>();
export function cssColor(root: HTMLElement, name: string, fallback: string): string {
  const cached = cssColorCache.get(name);
  if (cached !== undefined)
    return cached;
  const resolved = getComputedStyle(root).getPropertyValue(name).trim() || fallback;
  cssColorCache.set(name, resolved);
  return resolved;
}

/** Palette entries are keyed by name only, so a theme switch has to drop the whole cache. */
export function clearCssColorCache(): void {
  cssColorCache.clear();
}

export function argbToRgba(argb: number): [number, number, number, number] {
  return [DG.Color.r(argb) / 255, DG.Color.g(argb) / 255, DG.Color.b(argb) / 255, DG.Color.a(argb) / 255];
}

/** Draw a molecule with the ARGB background baked into the RDKit draw call — a colour applied
 *  afterwards would leave a pale anti-aliasing fringe around every bond. */
export function paintMoleculeOnColor(canvas: HTMLCanvasElement, smiles: string, w: number, h: number,
  argb: number): void {
  try {
    drawMoleculeToCanvas(0, 0, w, h, canvas, smiles, null,
      {normalizeDepiction: true, straightenDepiction: true}, null,
      {clearBackground: true, backgroundColour: argbToRgba(argb)});
  } catch (e) {
    // A structure RDKit cannot draw leaves the canvas at its background rather than failing the whole
    // pane — one bad row must not blank a navigator full of good ones.
  }
}

export function renderMoleculeOnColor(smiles: string, w: number, h: number, argb: number): HTMLElement {
  const canvas = ui.canvas(w, h);
  if (smiles)
    paintMoleculeOnColor(canvas as HTMLCanvasElement, smiles, w, h, argb);
  return canvas;
}

/** A filter panel over a hidden frame: the platform supplies the widgets, and the frame's own filter is
 *  read back as the set of keys that survive. Built on first open — an unopened filter would otherwise
 *  cost a fingerprinted column per matrix. */
export class FrameFilter<K> {
  private frame: DG.DataFrame | null = null;
  private view: DG.TableView | null = null;
  private group: DG.FilterGroup | null = null;
  private sub: Subscription | null = null;
  private keys: K[] = [];
  private pass: Set<K> | null = null;

  /** `configure` adds the filters the default set omits; `onChange` runs after every filter change. */
  constructor(
    private readonly build: () => {frame: DG.DataFrame, keys: K[]},
    private readonly configure: ((group: DG.FilterGroup) => void) | null,
    private readonly onChange: () => void) {}

  get active(): boolean {
    return this.pass !== null;
  }

  passes(key: K): boolean {
    return this.pass === null || this.pass.has(key);
  }

  root(): HTMLElement {
    if (this.group === null) {
      const built = this.build();
      this.frame = built.frame;
      this.keys = built.keys;
      this.view = DG.TableView.create(this.frame, false);
      this.group = this.view.getFiltersGroup();
      this.configure?.(this.group);
      this.sub = DG.debounce(this.frame.onFilterChanged, 300).subscribe(() => this.sync());
    }
    return this.group.root;
  }

  private sync(): void {
    const frame = this.frame;
    if (frame === null)
      return;
    if (frame.filter.trueCount === frame.rowCount)
      this.pass = null;
    else {
      const pass = new Set<K>();
      for (let i = 0; i < frame.rowCount; i++) {
        if (frame.filter.get(i))
          pass.add(this.keys[i]);
      }
      this.pass = pass;
    }
    this.onChange();
  }

  /** Rebuild against current data; stale rows key the previous set. The view is detached, not just
   *  dropped — the reference alone leaks the Dart-backed view. */
  reset(): void {
    this.sub?.unsubscribe();
    this.sub = null;
    this.view?.detach();
    this.view = null;
    this.frame = null;
    this.group = null;
    this.keys = [];
    this.pass = null;
  }
}
