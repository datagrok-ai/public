import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

/** Custom event id for toggling monomer highlights in macromolecule grid cells (HELM renderer
 *  and generic sequence renderer). Payload: {@link MacromoleculeHighlightEventArgs}. Fire via
 *  {@link grok.events.fireCustomEvent} or the convenience helpers in this module. */
export const MACROMOLECULE_HIGHLIGHT_EVENT_ID = 'bio-macromolecule-monomer-highlight';

/** Slot name on `column.temp` where per-row monomer highlights are stored by the macromolecule
 *  cell renderers. Shape: `Map<rowIdx, MacromoleculeHighlightEntry>`. */
export const MACROMOLECULE_HIGHLIGHT_TEMP = 'bio.macromolecule.highlightMonomers';

/** Forced alpha (0..1) applied to every highlight color — solid colors are not allowed. */
export const MACROMOLECULE_HIGHLIGHT_ALPHA = 0.4;

/** Default fill color (translucent yellow) as a DG.Color ARGB int. The alpha component of the
 *  input is ignored by the renderer and replaced with {@link MACROMOLECULE_HIGHLIGHT_ALPHA}. */
export const DEFAULT_MACROMOLECULE_HIGHLIGHT_FILL = DG.Color.argb(255, 255, 220, 0);
/** Default stroke color (amber) as a DG.Color ARGB int. */
export const DEFAULT_MACROMOLECULE_HIGHLIGHT_STROKE = DG.Color.argb(255, 230, 140, 0);

/** A single highlight record stored in `column.temp[MACROMOLECULE_HIGHLIGHT_TEMP]`, keyed by
 *  row index. */
export type MacromoleculeHighlightEntry = {
  /** 0-based monomer indices to highlight for this row.
   *  - HELM renderer: indices into `editor.m.atoms` of the rendered molecule.
   *  - Generic sequence renderer: 0-based position in the split sequence. */
  monomers: number[];
  /** DG.Color ARGB integer for the marker fill. Alpha is always replaced with
   *  {@link MACROMOLECULE_HIGHLIGHT_ALPHA}. Falls back to {@link DEFAULT_MACROMOLECULE_HIGHLIGHT_FILL}. */
  fillColor?: number;
  /** DG.Color ARGB integer for the marker outline (HELM renderer ring stroke). Alpha is always
   *  replaced with {@link MACROMOLECULE_HIGHLIGHT_ALPHA}. Falls back to
   *  {@link DEFAULT_MACROMOLECULE_HIGHLIGHT_STROKE}. Ignored by background-painting renderers. */
  strokeColor?: number;
};

export type MacromoleculeHighlightEventArgs = {
  /** Match by dataframe name. Optional — if omitted, `tableId` or no filter is used. */
  tableName?: string;
  /** Match by dataframe id. Optional. */
  tableId?: string;
  /** Name of the macromolecule column. */
  columnName: string;
  /** Row index in the dataframe (as stored, i.e. `gridCell.tableRowIndex`). */
  rowIdx: number;
  /** 0-based monomer indices to highlight. Pass `null` or `[]` to clear the row. */
  monomers: number[] | null;
  /** Optional DG.Color ARGB int for the marker fill. Alpha is forced to
   *  {@link MACROMOLECULE_HIGHLIGHT_ALPHA} to prevent solid colors. Ignored when clearing. */
  fillColor?: number;
  /** Optional DG.Color ARGB int for the marker outline (HELM ring). Alpha is forced to
   *  {@link MACROMOLECULE_HIGHLIGHT_ALPHA}. Ignored when clearing. */
  strokeColor?: number;
};

/** Optional color overrides accepted by the utility helpers (DG.Color ARGB ints). */
export type MacromoleculeHighlightColors = {
  fillColor?: number;
  strokeColor?: number;
};

/** Converts a DG.Color ARGB int to an `rgba(r,g,b,alpha)` CSS string. The input's alpha is
 *  discarded — the fixed {@link MACROMOLECULE_HIGHLIGHT_ALPHA} is used unless an explicit
 *  override is passed (renderers may use a smaller alpha for backgrounds). */
export function macromoleculeHighlightColorToCss(color: number, alpha?: number): string {
  return `rgba(${DG.Color.r(color)},${DG.Color.g(color)},${DG.Color.b(color)},${alpha ?? MACROMOLECULE_HIGHLIGHT_ALPHA})`;
}

/** Common input for the highlight helpers: accepts either a dataframe instance or its
 *  identifying name/id, plus either a column instance or its name. */
export type MacromoleculeHighlightTarget = {
  table?: DG.DataFrame;
  tableName?: string;
  tableId?: string;
  column?: DG.Column;
  columnName?: string;
};

function resolveTarget(t: MacromoleculeHighlightTarget): {tableName?: string, tableId?: string, columnName: string} {
  const columnName = t.columnName ?? t.column?.name;
  if (!columnName)
    throw new Error('MacromoleculeHighlight: columnName (or column) must be provided.');
  const tableName = t.tableName ?? t.table?.name;
  const tableId = t.tableId ?? t.table?.id;
  return {tableName, tableId, columnName};
}

/** Fires the macromolecule highlight event. Updates on the renderer side are applied
 *  asynchronously via the platform's custom-event bus; visible cells in the target column will
 *  be repainted. */
export function fireMacromoleculeHighlight(args: MacromoleculeHighlightEventArgs): void {
  grok.events.fireCustomEvent(MACROMOLECULE_HIGHLIGHT_EVENT_ID, args);
}

/** Highlight the given monomers in a single row. Replaces any previous highlight for that row.
 *  Optional `colors` override the defaults. */
export function setMacromoleculeMonomerHighlight(
  target: MacromoleculeHighlightTarget, rowIdx: number, monomers: number[],
  colors?: MacromoleculeHighlightColors,
): void {
  const {tableName, tableId, columnName} = resolveTarget(target);
  fireMacromoleculeHighlight({
    tableName, tableId, columnName, rowIdx, monomers,
    fillColor: colors?.fillColor, strokeColor: colors?.strokeColor,
  });
}

/** Clear the highlight for a single row. */
export function clearMacromoleculeMonomerHighlight(
  target: MacromoleculeHighlightTarget, rowIdx: number,
): void {
  const {tableName, tableId, columnName} = resolveTarget(target);
  fireMacromoleculeHighlight({tableName, tableId, columnName, rowIdx, monomers: null});
}

/** Per-row entry for {@link setMacromoleculeMonomerHighlights}: either `null` (clear), a plain
 *  monomer array (use defaults), or an object bundling monomers with optional color overrides. */
export type MacromoleculeHighlightRowSpec =
  | null
  | number[]
  | ({monomers: number[]} & MacromoleculeHighlightColors);

/** Set highlights for many rows at once. Each value is either `null` (clear that row), an
 *  array of monomer indices (use default colors), or `{monomers, fillColor?, strokeColor?}`. */
export function setMacromoleculeMonomerHighlights(
  target: MacromoleculeHighlightTarget,
  rowHighlights:
    | Map<number, MacromoleculeHighlightRowSpec>
    | Record<number, MacromoleculeHighlightRowSpec>,
): void {
  const {tableName, tableId, columnName} = resolveTarget(target);
  const entries: Array<[number, MacromoleculeHighlightRowSpec]> = rowHighlights instanceof Map ?
    Array.from(rowHighlights.entries()) :
    Object.entries(rowHighlights).map(([k, v]) => [Number(k), v] as [number, MacromoleculeHighlightRowSpec]);
  for (const [rowIdx, spec] of entries) {
    if (spec == null) {
      fireMacromoleculeHighlight({tableName, tableId, columnName, rowIdx, monomers: null});
    } else if (Array.isArray(spec)) {
      fireMacromoleculeHighlight({tableName, tableId, columnName, rowIdx, monomers: spec});
    } else {
      fireMacromoleculeHighlight({
        tableName, tableId, columnName, rowIdx,
        monomers: spec.monomers, fillColor: spec.fillColor, strokeColor: spec.strokeColor,
      });
    }
  }
}

/** Clears all highlights stored on the column by re-firing a clear event for every currently-
 *  highlighted row. */
export function clearAllMacromoleculeMonomerHighlights(target: MacromoleculeHighlightTarget): void {
  const {tableName, tableId, columnName} = resolveTarget(target);
  const col = target.column ?? target.table?.columns.byName(columnName);
  if (!col) return;
  const map = col.temp[MACROMOLECULE_HIGHLIGHT_TEMP] as Map<number, MacromoleculeHighlightEntry> | undefined;
  if (!map || map.size === 0) return;
  const rows = Array.from(map.keys());
  for (const rowIdx of rows)
    fireMacromoleculeHighlight({tableName, tableId, columnName, rowIdx, monomers: null});
}
