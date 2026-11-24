import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

// The platform keeps at most 20 MB per undo step: this many cells of kept or copied columns stay under it.
const MAX_KEPT_CELLS = 2_500_000;

function push(name: string, context: DG.DataFrame | DG.ViewBase, undo: () => void, redo?: () => void): boolean {
  DG.UndoService?.push(name, undo, {redo, context});
  return DG.UndoService != null;
}

/** Table changes that Ctrl+Z takes back, for the AI's tools and scripts: each makes the change and records the
 * step. A boolean result is false when the change was made without a step, so it cannot be undone. */
export const undoable = {
  removeColumns(t: DG.DataFrame, names: string[]): boolean {
    const all = t.columns.names();
    const removed = [...new Set(names.map((n) => t.getCol(n)))].map((col) => ({col, at: all.indexOf(col.name)}))
      .sort((a, b) => a.at - b.at);
    const remove = () => removed.forEach(({col}) => t.columns.remove(col.name));
    remove();
    if (removed.length === 0 || removed.length * t.rowCount > MAX_KEPT_CELLS)
      return false;
    return push(removed.length === 1 ? `Remove column ${removed[0].col.name}` : `Remove ${removed.length} columns`, t,
      () => removed.forEach(({col, at}) => t.columns.insert(col, at)), remove);
  },

  /** Removes the rows [where] holds for with the platform's command, which records the step itself and drops it
   * when the rows are too large to keep. The user's selection of the other rows stays. */
  async removeRows(t: DG.DataFrame, where: (i: number) => boolean): Promise<boolean> {
    if (grok.shell.t !== t)
      throw new Error(`Rows are removed from the current table only: make the view of "${t.name}" current first`);
    const rowCount = t.rowCount;
    const mask = Array.from({length: rowCount}, (_, i) => !!where(i));
    const removed = mask.filter(Boolean).length;
    if (removed === 0)
      return false;
    const selected = t.selection.clone();
    t.selection.init((i) => mask[i]);
    try {
      await grok.functions.call('CmdRemoveSelectedRows');
    } finally {
      if (t.rowCount === rowCount - removed) {
        const kept = Array.from({length: rowCount}, (_, i) => i).filter((i) => !mask[i]);
        t.selection.init((i) => selected.get(kept[i]));
      } else
        t.selection.copyFrom(selected);
    }
    return DG.UndoService?.undoName === `Remove ${removed} rows`;
  },

  /** Records a column the script has just added to its table. */
  columnAdded(col: DG.Column): boolean {
    const t = col.dataFrame;
    const at = t.columns.names().indexOf(col.name);
    return push(`Add column ${col.name}`, t, () => t.columns.remove(col.name), () => t.columns.insert(col, at));
  },

  /** Runs [change], which rewrites the values of [col]; a column too large to copy is changed without a step.
   * A change that throws is rolled back. */
  async changeValues(col: DG.Column, change: () => any): Promise<boolean> {
    const before = col.length <= MAX_KEPT_CELLS ? col.clone() : null;
    const restore = () => col.dataFrame.columns.replace(col.name, before!);
    try {
      await change();
    } catch (e) {
      if (before)
        restore();
      throw e;
    }
    return before != null && push(`Change ${col.name}`, col.dataFrame, restore);
  },

  /** Runs [change], which renames [col] or sets its tags: format, color coding, units, semantic type.
   * A change that throws is rolled back. */
  async editColumn(col: DG.Column, change: () => any): Promise<boolean> {
    const state = () => ({name: col.name, tags: Object.fromEntries(col.tags) as Record<string, string>});
    const restore = (to: ReturnType<typeof state>) => {
      col.name = to.name;
      for (const key of Object.keys(state().tags)) {
        if (!(key in to.tags))
          delete col.tags[key];
      }
      Object.assign(col.tags, to.tags);
    };
    const before = state();
    try {
      await change();
    } catch (e) {
      restore(before);
      throw e;
    }
    const after = state();
    return JSON.stringify(before) !== JSON.stringify(after) &&
      push(`Edit column ${after.name}`, col.dataFrame, () => restore(before), () => restore(after));
  },

  addViewer(view: DG.TableView, type: string, options?: object): DG.Viewer {
    let viewer = view.addViewer(type, options);
    push(`Add ${type.toLowerCase()}`, view, () => {
      if (Array.from(view.viewers).includes(viewer))
        viewer.close();
    }, () => viewer = view.addViewer(type, options));
    return viewer;
  },
};
