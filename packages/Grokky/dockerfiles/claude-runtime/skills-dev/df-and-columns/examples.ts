/**
 * Compile-time verification of every fenced example in SKILL.md.
 *
 * Every `datagrok-exec` block referenced from SKILL.md appears here as a
 * named function with the same code. If any example uses a non-existent
 * method or wrong arg shape, this file fails to type-check.
 *
 * Each example is wrapped in a function with the same global signature the
 * runtime injects into a `datagrok-exec` block: `(grok, ui, DG, view, t)`.
 * The runtime does NOT inject a `grokky` namespace — these examples must
 * compile against raw js-api only.
 */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

type ExampleFn = (
  grok: typeof import('datagrok-api/grok'),
  ui: typeof import('datagrok-api/ui'),
  DG: typeof import('datagrok-api/dg'),
  view: DG.ViewBase,
  t: DG.DataFrame,
) => Promise<unknown> | unknown;

export const examples: Array<{name: string; fn: ExampleFn}> = [
  // ---------- Finding columns ----------

  {name: 'find-molecule-column', fn: async (_grok, ui, DG, _view, t) => {
    const molCol = t.columns.bySemType(DG.SEMTYPE.MOLECULE);
    return ui.divText(molCol ? `Molecule column: ${molCol.name}` : 'No molecule column');
  }},

  {name: 'find-all-molecule-columns', fn: async (_grok, ui, DG, _view, t) => {
    const molCols = t.columns.bySemTypeAll(DG.SEMTYPE.MOLECULE);
    return ui.divV(molCols.map((c) => ui.divText(c.name)));
  }},

  {name: 'find-mw-fallback-substring', fn: async (_grok, ui, _DG, _view, t) => {
    const wanted = 'mw';
    let col: DG.Column | null | undefined = t.col(wanted);
    if (!col) {
      const lc = wanted.toLowerCase();
      col = t.columns.firstWhere((c) => c.name.toLowerCase().includes(lc));
    }
    return ui.divText(col ? col.name : 'no match');
  }},

  // ---------- Describing a column ----------

  {name: 'describe-one-column', fn: async (_grok, ui, _DG, _view, t) => {
    const col = t.getCol('age');
    const s = col.stats;
    return ui.tableFromMap({
      name: col.name,
      type: col.type,
      semType: col.semType ?? '(none)',
      length: col.length,
      missing: s.missingValueCount,
      unique: s.uniqueCount,
      min: s.min, max: s.max, avg: s.avg, stdev: s.stdev,
    });
  }},

  {name: 'describe-all-numeric', fn: async (_grok, _ui, DG, _view, t) => {
    const rows = Array.from(t.columns.numerical, (c) => ({
      name: c.name, min: c.stats.min, max: c.stats.max, avg: c.stats.avg,
    }));
    return DG.Viewer.grid(DG.DataFrame.fromObjects(rows)!).root;
  }},

  {name: 'top-n-categories', fn: async (_grok, ui, _DG, _view, t) => {
    const col = t.getCol('class');
    const counts = new Map<string, number>();
    for (let i = 0; i < col.length; i++) {
      const v = col.get(i) as string | null;
      if (v === null || v === undefined) continue;
      counts.set(v, (counts.get(v) ?? 0) + 1);
    }
    const top = Array.from(counts.entries())
      .sort((a, b) => b[1] - a[1])
      .slice(0, 5);
    return ui.tableFromMap(Object.fromEntries(top));
  }},

  // ---------- Cloning ----------

  {name: 'clone-full', fn: async (_grok, _ui, DG, _view, t) => {
    const copy = t.clone();
    return DG.Viewer.grid(copy).root;
  }},

  {name: 'clone-filtered-subset', fn: async (_grok, _ui, DG, _view, t) => {
    const small = t.clone(t.filter, ['smiles', 'activity']);
    return DG.Viewer.grid(small).root;
  }},

  // ---------- Adding columns ----------

  {name: 'add-empty-float-with-meta', fn: async (_grok, ui, DG, _view, t) => {
    const col = t.columns.addNewFloat('Ki');
    col.semType = DG.SEMTYPE.Ki;
    col.meta.units = 'nM';
    col.meta.format = '0.00';
    return ui.divText(`Added ${col.name} (${col.type})`);
  }},

  {name: 'add-from-values', fn: async (_grok, ui, DG, _view, t) => {
    const values = [1.2, 3.4, null, 5.6];
    const col = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'score', values);
    t.columns.add(col);
    return ui.divText(`Added ${col.name}`);
  }},

  {name: 'add-virtual-column', fn: async (_grok, ui, DG, _view, t) => {
    const labels = t.columns.addNewVirtual('Label',
      (i: number) => `${t.get('compound', i)}@${t.get('target', i)}`,
      DG.TYPE.STRING);
    return ui.divText(`Added virtual column ${labels.name}`);
  }},

  {name: 'add-init-by-index', fn: async (_grok, ui, _DG, _view, t) => {
    const col = t.columns.addNewInt('rowSquared');
    col.init((i: number) => i * i);
    return ui.divText(`Initialized ${col.name}`);
  }},

  {name: 'add-with-unused-name', fn: async (_grok, ui, _DG, _view, t) => {
    const name = t.columns.getUnusedName('Ki');
    const col = t.columns.addNewFloat(name);
    return ui.divText(`Final name: ${col.name}`);
  }},

  // ---------- Removing & renaming ----------

  {name: 'remove-one-column', fn: async (_grok, _ui, _DG, _view, t) => {
    t.columns.remove('temp');
  }},

  {name: 'remove-bulk-loop', fn: async (_grok, _ui, _DG, _view, t) => {
    for (const name of ['_tmp1', '_tmp2', 'temp3'])
      if (t.columns.contains(name))
        t.columns.remove(name);
  }},

  {name: 'remove-calculated', fn: async (_grok, ui, _DG, _view, t) => {
    const names = t.columns.toList().filter((c) => c.meta.formula != null).map((c) => c.name);
    for (const n of names) t.columns.remove(n);
    return ui.divText(`Removed: ${names.join(', ') || '(none)'}`);
  }},

  {name: 'rename-column', fn: async (_grok, _ui, _DG, _view, t) => {
    t.getCol('pIC50').name = 'pIC50_old';
  }},

  // ---------- Metadata ----------

  {name: 'set-semtype', fn: async (_grok, _ui, DG, _view, t) => {
    t.getCol('smiles').semType = DG.SEMTYPE.MOLECULE;
  }},

  {name: 'set-friendly-name', fn: async (_grok, _ui, _DG, _view, t) => {
    t.getCol('pIC50').meta.friendlyName = 'Potency';
  }},

  {name: 'set-meta-bulk', fn: async (_grok, _ui, DG, _view, t) => {
    const col = t.getCol('Ki');
    col.semType = DG.SEMTYPE.Ki;
    col.meta.units = 'nM';
    col.meta.format = '0.00';
    col.meta.description = null;
  }},

  {name: 'clear-all-user-meta', fn: async (_grok, _ui, _DG, _view, t) => {
    const col = t.getCol('activity');
    col.meta.friendlyName = null;
    col.meta.units = null;
    col.meta.format = null;
    col.meta.description = null;
  }},

  // ---------- Color coding ----------

  {name: 'color-coding-linear', fn: async (_grok, _ui, _DG, _view, t) => {
    t.getCol('age').meta.colors.setLinear(
      ['#ff0000', '#ffff00', '#00ff00'] as unknown as number[],
      {min: 19, max: 70},
    );
  }},

  {name: 'color-coding-categorical', fn: async (_grok, _ui, _DG, _view, t) => {
    t.getCol('race').meta.colors.setCategorical(
      {'Asian': '#0000FF', 'Black': '#FF0000'},
      {fallbackColor: '#CCCCCC'},
    );
  }},

  {name: 'color-coding-conditional', fn: async (_grok, _ui, _DG, _view, t) => {
    t.getCol('height').meta.colors.setConditional({
      '20-170': '#00FF00',
      '170-190': '#220505',
    });
  }},

  {name: 'color-coding-off', fn: async (_grok, _ui, _DG, _view, t) => {
    t.getCol('activity').meta.colors.setDisabled();
  }},

  // ---------- Cell values ----------

  {name: 'bulk-init-via-fn', fn: async (_grok, _ui, _DG, _view, t) => {
    const col = t.columns.addNewInt('rowSquared');
    col.init((i: number) => i * i);
  }},

  {name: 'hot-loop-with-fire', fn: async (_grok, _ui, _DG, _view, t) => {
    const src = t.getCol('input');
    const dst = t.getCol('output');
    for (let i = 0; i < t.rowCount; i++)
      dst.set(i, (src.get(i) as number) * 2, false);
    dst.fireValuesChanged();
  }},

  {name: 'correct-typed-column-via-constant', fn: async (_grok, _ui, DG, _view, t) => {
    t.columns.addNew('x', DG.COLUMN_TYPE.FLOAT);
    t.columns.addNewFloat('y');
  }},

  // Silence "imported but never used" while keeping the imports for type checking.
  {name: '_noop-imports', fn: () => {
    void grok; void ui; void DG;
  }},
];
