import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

/** What the platform's own pickers offer, read one way by every u2/dg picker: the open tables a picker can name, and
 * the columns of a frame a picker lists. */
export class Platform {
  /** The names of the open tables, in the shell's order and each once. A table renamed to null is listed by the
   * shell as `undefined` (SEAM S-23) and one named '' shows nothing, so neither is offered. */
  static tableNames(): string[] {
    return [...new Set(grok.shell.tableNames.filter((name) => Platform._named(name)))];
  }

  /** Whether `a` and `b` are one frame: js-api wraps a frame anew on every read, so the dart handle (or, failing
   * that, the name) is its identity. */
  static sameTable(a: DG.DataFrame | null, b: DG.DataFrame | null): boolean {
    if (a === b)
      return true;
    if (a == null || b == null)
      return false;
    return a.dart != null && b.dart != null ? a.dart === b.dart : a.name === b.name;
  }

  /** Whether `table` is open in the workspace. */
  static isOpen(table: DG.DataFrame): boolean {
    return (grok.shell.tables ?? []).some((open) => Platform.sameTable(open, table));
  }

  /** The columns of `table` the platform's pickers list (`column_grid.dart:909-911`), in table order: none hidden
   * behind `~` (Chem's `~<column>.canonicalSmiles` cache) and none without a name; `filter` narrows them further. */
  static listedColumns(table: DG.DataFrame, filter?: (column: DG.Column) => boolean): DG.Column[] {
    return table.columns.toList().filter((c) => Platform._named(c.name) && !c.name.startsWith('~') &&
      (filter === undefined || filter(c)));
  }

  private static _named(name: unknown): boolean {
    return typeof name === 'string' && name !== '';
  }
}
