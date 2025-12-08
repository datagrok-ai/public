/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {DBValueObject, EntryPointOptions, SchemaAndConnection} from './types';
import {DBExplorerRenderer} from './renderer';


export class DBExplorerObjectHandler extends DG.ObjectHandler {
  get type(): string {
    return 'db-explorer-value';
  }
  private renderer: DBExplorerRenderer;
  isApplicable(x: any): boolean {
    return x instanceof DBValueObject;
  }

  constructor(
    public options: EntryPointOptions,
    private schemaInfoPromise: () => Promise<SchemaAndConnection | null>,
    schemaName: string
  ) {
    super();
    this.renderer = new DBExplorerRenderer(schemaInfoPromise, schemaName, options);
  }

  addCustomRenderer(
    check: (tableName: string, colName: string, value: string | number) => boolean,
    renderer: (value: string | number, connection: DG.DataConnection) => HTMLElement
  ) {
    this.renderer.addCustomRenderer(check, renderer);
  }

  addDefaultHeaderReplacerColumns(columns: string[]) {
    this.renderer.addDefaultHeaderReplacerColumns(columns);
    return this;
  }

  addCustomSelectedColumns(columns: {[tableName: string]: string[]}) {
    this.renderer.addCustomSelectedColumns(columns);
  }

  addUniqueColumns(columns: {[tableName: string]: string}) {
    this.renderer.addUniqueColNames(columns);
  }
  addHeaderReplacers(replacers: {[tableName: string]: string}) {
    this.renderer.addHeaderReplacers(replacers);
  }

  renderInnerCard(x: DBValueObject) {
    return ui.wait(async () => {
      return await this.renderer.renderTable(
        x.table,
        x.column,
        this.options.valueConverter(x.value ?? ''),
        this.options.joinOptions
      );
    });
  }

  renderCard(x: DBValueObject, _context?: any): HTMLElement {
    const c = ui.card(this.renderInnerCard(x));
    c.style.width = 'unset';
    c.addEventListener('click', () => (grok.shell.o = x));
    return c;
  }

  renderTooltip(x: DBValueObject, _context?: any): HTMLElement {
    return this.renderInnerCard(x);
  }

  renderInnerProperties(tableName: string, x: DBValueObject, paneName: string) {
    const acc = ui.accordion(tableName);
    acc.addPane(paneName, () => this.renderInnerCard(x), true);
    this.renderer
      .getTable(tableName, x.column, this.options.valueConverter(x.value ?? ''), this.options.joinOptions)
      .then((df) => {
        this.renderer.renderAssociations(acc, this.schemaInfoPromise, tableName, df);
      });
    return acc;
  }

  renderProperties(x: DBValueObject, _context?: any): HTMLElement {
    const acc = this.renderInnerProperties(x.table, x, this.options.valueConverter(x.value ?? '').toString());
    if (x.semValue) {
      const origAcc = ui.panels.infoPanel(x.semValue);
      origAcc.context = x.semValue;
      if (x.semValue.semType)
        x.semValue.tags[DG.Tags.Quality] = x.semValue.semType;
      origAcc?.end();

      return ui.divV([acc.root, origAcc.root]);
    }
    return acc.root;
  }
}

export class SemValueObjectHandler extends DBExplorerObjectHandler {
  get type(): string {
    return this.semanticType;
  }

  private _rgExample: (typeof DG.ObjectHandler.prototype.regexpExample) | null = null;

  override get regexpExample() {
    return this._rgExample;
  }
  isApplicable(x: any): boolean {
    return x instanceof DG.SemanticValue && x.semType == this.semanticType;
  }

  get entryPoint() {
    return {table: this.tableName, column: this.columnName};
  }

  constructor(
    private semanticType: string,
    private tableName: string,
    private columnName: string,
    options: EntryPointOptions,
    schemaInfoPromise: () => Promise<SchemaAndConnection | null>,
    schemaName: string
  ) {
    super(options, schemaInfoPromise, schemaName);
    if (options && options.regexpExample && options.regexpExample.nonVariablePart && options.regexpExample.regexpMarkup)
      this._rgExample = options.regexpExample;
  }

  renderCard(x: any, context?: any): HTMLElement {
    return super.renderCard(new DBValueObject(this.tableName, this.columnName, this.options.valueConverter(x.value ?? ''), x), context);
  }

  renderTooltip(x: any, context?: any): HTMLElement {
    return super.renderTooltip(new DBValueObject(this.tableName, this.columnName, this.options.valueConverter(x.value ?? ''), x), context);
  }

  renderProperties(x: any, context?: any): HTMLElement {
    return super.renderProperties(new DBValueObject(this.tableName, this.columnName, this.options.valueConverter(x.value ?? ''), x), context);
  }
}
