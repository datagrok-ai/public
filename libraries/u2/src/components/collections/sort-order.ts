/* A table's order as one string, the way DataTable, BasicTable and CatalogTable keep it: 'name' ascending, '!name'
   descending, '' none; keys after a comma break ties ('name,!date'). */

export interface SortKey {
  column: string;
  descending: boolean;
}

export class SortOrder {
  /** The first key of `sort`; null for none. */
  static parse(sort: string): SortKey | null {
    const first = sort.split(',')[0];
    const descending = first.startsWith('!');
    const column = descending ? first.slice(1) : first;
    return column === '' ? null : {column, descending};
  }

  /** `column` ascending, or descending. */
  static of(column: string, descending = false): string {
    return descending ? `!${column}` : column;
  }

  /** What a press on `column`'s header makes of `sort`: ascending, then descending, then none. */
  static next(sort: string, column: string): string {
    const key = SortOrder.parse(sort);
    return key?.column !== column ? column : key.descending ? '' : `!${column}`;
  }
}
