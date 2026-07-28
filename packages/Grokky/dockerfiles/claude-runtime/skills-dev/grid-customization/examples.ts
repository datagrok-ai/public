/**
 * Compile-time verification of every fenced example in SKILL.md.
 *
 * Each example is a named function with the same signature the runtime
 * injects into a `datagrok-exec` block: `(grok, ui, DG, view, t)`. No
 * `grokky` global exists — examples call the raw js-api directly.
 *
 * If any example uses a non-existent method or the wrong arg shape, this
 * file fails to type-check. No `@ts-ignore` / `as any` cheating: where the
 * js-api type is narrower than the platform runtime actually accepts (e.g.
 * `setLinear` is typed `number[]` but accepts hex strings too — see
 * `ApiSamples/scripts/grid/color-coding/color-coding.js`), the example
 * uses the typed shape. The runtime accepts strings; the cast is loud and
 * one-line, not a blanket escape hatch.
 */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

type ExampleFn = (
  grok: typeof import('datagrok-api/grok'),
  ui: typeof import('datagrok-api/ui'),
  DG: typeof import('datagrok-api/dg'),
  view: DG.TableView,
  t: DG.DataFrame,
) => Promise<unknown> | unknown;

// Hex-string palettes are accepted by setLinear at runtime even though the
// public type signature is `number[]`. Single localized cast — not a global
// loosening.
const asColorList = (xs: string[]): number[] => xs as unknown as number[];

export const examples: Array<{name: string; fn: ExampleFn}> = [

  // ---------------------------------------------------------------------
  // Sorting
  // ---------------------------------------------------------------------

  {name: 'sort-activity-desc', fn: (_g, _u, _DG, view, _t) => {
    view.grid.sort(['activity'], [false]);
  }},

  {name: 'sort-class-then-activity', fn: (_g, _u, _DG, view, _t) => {
    view.grid.sort(['class', 'activity'], [true, false]);
  }},

  {name: 'clear-sort', fn: (_g, _u, _DG, view, _t) => {
    view.grid.sort([], []);
  }},

  // ---------------------------------------------------------------------
  // Visibility & order
  // ---------------------------------------------------------------------

  {name: 'hide-single-column', fn: (_g, _u, _DG, view, _t) => {
    const gc = view.grid.col('_internalId');
    if (gc) gc.visible = false;
  }},

  {name: 'show-only-five', fn: (_g, _u, _DG, view, _t) => {
    view.grid.columns.setVisible(['SMILES', 'MW', 'LogP', 'activity', 'class']);
  }},

  {name: 'reorder-front-three', fn: (_g, _u, _DG, view, _t) => {
    view.grid.columns.setOrder(['SMILES', 'MW', 'LogP']);
  }},

  {name: 'reveal-all-columns', fn: (_g, _u, _DG, view, _t) => {
    const cols = view.grid.columns;
    for (let i = 1; i < cols.length; i++) {
      const gc = cols.byIndex(i);
      if (gc) gc.visible = true;
    }
  }},

  // ---------------------------------------------------------------------
  // Widths
  // ---------------------------------------------------------------------

  {name: 'widths-by-name', fn: (_g, _u, _DG, view, _t) => {
    const smiles = view.grid.col('SMILES'); if (smiles) smiles.width = 250;
    const mw = view.grid.col('MW'); if (mw) mw.width = 80;
    const logp = view.grid.col('LogP'); if (logp) logp.width = 70;
  }},

  {name: 'width-policy-optimal', fn: (_g, _u, DG, view, _t) => {
    view.grid.setColumnsWidthType(DG.ColumnWidthType.Optimal);
  }},

  // ---------------------------------------------------------------------
  // Pinning & freezing
  // ---------------------------------------------------------------------

  {name: 'pin-smiles', fn: (_g, _u, _DG, view, _t) => {
    const gc = view.grid.col('SMILES');
    if (gc) gc.pin();
  }},

  {name: 'pin-two', fn: (_g, _u, _DG, view, _t) => {
    const smiles = view.grid.col('SMILES'); if (smiles) smiles.pin();
    const id = view.grid.col('ID'); if (id) id.pin();
  }},

  {name: 'unpin-smiles', fn: (_g, _u, _DG, view, _t) => {
    const gc = view.grid.col('SMILES');
    if (gc) gc.unpin();
  }},

  {name: 'freeze-first-two', fn: (_g, _u, _DG, view, _t) => {
    view.grid.setOptions({frozenColumns: 2});
  }},

  // ---------------------------------------------------------------------
  // Row height
  // ---------------------------------------------------------------------

  {name: 'row-height-100', fn: (_g, _u, _DG, view, _t) => {
    view.grid.setOptions({rowHeight: 100});
  }},

  // ---------------------------------------------------------------------
  // Formatting
  // ---------------------------------------------------------------------

  {name: 'format-ic50-and-count', fn: (_g, _u, _DG, _view, t) => {
    t.col('IC50')!.meta.format = '0.00';
    t.col('count')!.meta.format = '#,##0';
  }},

  // ---------------------------------------------------------------------
  // Color coding
  // ---------------------------------------------------------------------

  {name: 'color-linear-activity-green-red', fn: (_g, _u, _DG, _view, t) => {
    t.col('activity')!.meta.colors.setLinear(
      asColorList(['#00FF00', '#FFFF00', '#FF0000']),
    );
  }},

  {name: 'color-linear-min-max', fn: (_g, _u, _DG, _view, t) => {
    t.col('IC50')!.meta.colors.setLinear(
      asColorList(['#FF0000', '#FFFF00', '#00FF00']),
      {min: 0, max: 100, belowMinColor: '#000080', aboveMaxColor: '#000000'},
    );
  }},

  {name: 'color-categorical-class', fn: (_g, _u, _DG, _view, t) => {
    t.col('class')!.meta.colors.setCategorical(
      {A: '#FF0000', B: '#00FF00', C: '#0000FF'},
      {fallbackColor: '#CCCCCC'},
    );
  }},

  {name: 'color-conditional-height', fn: (_g, _u, _DG, _view, t) => {
    t.col('height')!.meta.colors.setConditional({
      '20-170': '#00FF00',
      '170-190': '#FFFF00',
      '190-': '#FF0000',
    });
  }},

  {name: 'color-off', fn: (_g, _u, _DG, _view, t) => {
    t.col('activity')!.meta.colors.setDisabled();
  }},

  // ---------------------------------------------------------------------
  // Reset (no single API — composed)
  // ---------------------------------------------------------------------

  {name: 'reset-grid-keep-viewers', fn: (_g, _u, DG, view, t) => {
    const cols = view.grid.columns;
    for (let i = 1; i < cols.length; i++) {
      const gc = cols.byIndex(i);
      if (gc) gc.visible = true;
    }
    view.grid.setColumnsWidthType(DG.ColumnWidthType.Optimal);
    view.grid.sort([], []);
    for (const c of t.columns.toList())
      c.meta.colors.setDisabled();
  }},

  // ---------------------------------------------------------------------
  // Combined HitTriage-style recipe
  // ---------------------------------------------------------------------

  {name: 'hit-list-multi-step', fn: (_g, _u, _DG, view, t) => {
    t.col('activity')!.meta.colors.setLinear(
      asColorList(['#FF0000', '#FFFF00', '#00FF00']),
    );
    const idx = view.grid.col('index'); if (idx) idx.visible = false;
    const smiles = view.grid.col('SMILES'); if (smiles) smiles.pin();
    view.grid.sort(['activity'], [false]);
  }},
];

// Reference the imported symbols so the bundle doesn't complain about
// unused imports — these globals are always present in datagrok-exec.
export const _refs = {grok, ui, DG};
