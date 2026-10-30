import type {IProperty, RowView, Scope} from '@datagrok-libraries/u2';
import {div, divV, span, text} from '@datagrok-libraries/u2';
import {Editors} from '@datagrok-libraries/u2/dg';
import type {DomainApp} from './app.js';
import {DomainTable} from './facade.js';
import type {DomainAppViewOptions} from './facade.js';

/** The pages `/domains/<schema>/<table>` opens for a table no app configured: the conditions beside the rows,
 * cards where a picture tells the rows apart and a table of columns where nothing does, and an entity page headed
 * by the row's name and its picture — its own, or the one a table extending it one to one holds (a compound's
 * structure kept beside it). */
export class GenericPage {
  static readonly HEADER_PICTURE = {width: 260, height: 170};
  /** The columns a table of rows starts with; Output adds the rest. */
  static readonly COLUMNS = 10;

  static options(handle: DomainTable): DomainAppViewOptions {
    const pictured = DomainTable.pictureOf(handle.properties) !== undefined;
    return {filters: {mode: 'builder', placement: 'left'}, setup: (app, self) => GenericPage._setup(app, self),
      ...pictured ? {mode: 'cards' as const, itemHeight: 96} : {view: 'table' as const}};
  }

  private static _setup(app: DomainApp, self: Scope): void {
    app.root.classList.add('u2-domain-generic');
    const columns = app.output.columns.peek();
    if (columns.length > GenericPage.COLUMNS)
      app.output.columns.value = columns.slice(0, GenericPage.COLUMNS);
    const picture = GenericPage._picture(app.table);
    self.effect(() => {
      const source = app.entitySource.value;
      const row = source?.currentRow.value ?? null;
      const page = app.root.querySelector('[data-u2-part="entity-page"]');
      if (row === null || source!.isDraft || page === null)
        return;
      const art = div([], 'u2-domain-header-picture');
      const key = app.table.properties.filter((p) => app.table.info.businessKey.includes(p.name!))
        .map((p) => `${DomainTable.captionOf(p)} ${text(row[p.name!])}`).join(' · ');
      const header = div([art, divV([span(app.table.renderer.caption(row), 'u2-domain-header-title'),
        span(key, 'u2-domain-header-key')])], 'u2-domain-header');
      header.dataset.u2Part = 'header';
      page.prepend(header);
      void picture.then((found) => found?.(row)).then((el) => el ? art.append(el) : art.remove());
      return () => header.remove();
    });
  }

  /** How a row's picture is drawn: from its own pictured column, else from the pictured table whose key is a link
   * to this one; null where neither holds one. */
  private static async _picture(table: DomainTable): Promise<((row: RowView) => Promise<HTMLElement | null>) | null> {
    const own = DomainTable.pictureOf(table.properties);
    if (own !== undefined)
      return async (row) => GenericPage._draw(own, row[own.name!]);
    const key = table.info.businessKey;
    if (key.length !== 1)
      return null;
    for (const child of table.info.childTables) {
      const extension = await DomainTable.shared(`${child.schema}.${child.table}`).catch(() => null);
      const prop = extension === null ? undefined : DomainTable.pictureOf(extension.properties);
      if (extension === null || prop === undefined || extension.info.businessKey.join() !== child.fkColumn)
        continue;
      return async (row) => {
        const filter = {property: child.fkColumn, operator: '=', value: String(row[key[0]])};
        const [found] = await extension.table.query({filter, limit: 1});
        return found === undefined ? null : GenericPage._draw(prop, found[prop.name!]);
      };
    }
    return null;
  }

  private static _draw(prop: IProperty, value: unknown): HTMLElement | null {
    return value === null || value === undefined || value === '' ? null :
      Editors.display(prop, value, GenericPage.HEADER_PICTURE);
  }
}
