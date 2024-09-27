/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {DBExplorerConfig, EntryPointOptions,
  QueryJoinOptions, ReferencedByObject, ReferenceObject, SchemaInfo} from './types';
import {DBExplorerObjectHandler, SemValueObjectHandler} from './object-handlers';

export class DBExplorer {
  private schemasLoaded: boolean = false;
  private connection: DG.DataConnection | null = null;
  private references: ReferenceObject = {};
  private referencedBy: ReferencedByObject = {};
  private _dbLoadPromise: Promise<void>;
  private objHandlers: DBExplorerObjectHandler[] = [];
  constructor(
    private connectionName: string,
    private schemaName: string
  ) {
    this._dbLoadPromise = this.loadDbSchema();
  }

  public get dbSchema(): Promise<SchemaInfo> {
    return this._dbLoadPromise.then(() => ({references: this.references, referencedBy: this.referencedBy}));
  }

  private async loadDbSchema() {
    this.connection = await grok.dapi.connections.filter(`name="${this.connectionName}"`).first();
    if (this.connection == null)
      throw new Error(`Connection ${this.connectionName} not found`);

    const tables = await grok.dapi.connections.getSchema(this.connection, this.schemaName);
    if (!tables)
      throw new Error(`Schema ${this.schemaName} not found`);

    tables.forEach((table) => {
      const tableName = table.friendlyName ?? table.name;
      this.references[tableName] = {};
      const t = this.references[tableName];
      table.columns.forEach((column) => {
        const ref = column.referenceInfo;
        if (ref && ref.table && ref.column) {
          t[column.name] = {refTable: ref.table, refColumn: ref.column};
          if (!this.referencedBy[ref.table])
            this.referencedBy[ref.table] = {};

          if (!this.referencedBy[ref.table][ref.column])
            this.referencedBy[ref.table][ref.column] = [];

          this.referencedBy[ref.table][ref.column].push({refTable: tableName, refColumn: column.name});
        }
      });
    });

    const handler = new DBExplorerObjectHandler(
      {valueConverter: (a) => a, joinOptions: []},
      {references: this.references, referencedBy: this.referencedBy},
      this.connection!,
      this.schemaName
    );
    DG.ObjectHandler.register(handler);
    this.objHandlers.push(handler);
    this.schemasLoaded = true;
  }

  public async addCustomRelation(tableName: string, columnName: string, refTable: string, refColumn: string) {
    if (!this.schemasLoaded)
      await this._dbLoadPromise;
    if (!this.references[tableName])
      this.references[tableName] = {};
    this.references[tableName][columnName] = {refTable, refColumn};
    if (!this.referencedBy[refTable])
      this.referencedBy[refTable] = {};
    this.referencedBy[refTable][refColumn] ??= [];
    this.referencedBy[refTable][refColumn].push({refTable: tableName, refColumn: columnName});
  }

  public async addEntryPoint(
    semanticType: string,
    tableName: string,
    columnName: string,
    options?: Partial<EntryPointOptions>
  ) {
    const schema = await this.dbSchema;
    const fullOpts = {...{valueConverter: (a: string | number) => a, joinOptions: []}, ...(options ?? {})};
    const handler = new SemValueObjectHandler(
      semanticType,
      tableName,
      columnName,
      fullOpts,
      schema,
      this.connection!,
      this.schemaName
    );
    DG.ObjectHandler.register(handler);
    this.objHandlers.push(handler);
    return this;
  }

  addJoinOptions(joinOpts: QueryJoinOptions[]) {
    this.objHandlers.forEach((handler) => handler.options.joinOptions.push(...joinOpts));
    return this;
  }

  public addCustomRenderer(
    check: (tableName: string, colName: string, value: string | number) => boolean,
    renderer: (value: string | number, connection: DG.DataConnection) => HTMLElement
  ) {
    this.objHandlers.forEach((handler) => handler.addCustomRenderer(check, renderer));
    return this;
  }

  public addHeaderReplacers(replacers: {[tableName: string]: string}) {
    this.objHandlers.forEach((handler) => handler.addHeaderReplacers(replacers));
    return this;
  }

  public addEntryPointValueConverter(func: (a: string | number) => string | number) {
    this.objHandlers.forEach((handler) => handler.options.valueConverter = func);
    return this;
  }

  public addDefaultHeaderReplacerColumns(columns: string[]) {
    this.objHandlers.forEach((handler) => handler.addDefaultHeaderReplacerColumns(columns));
    return this;
  }

  public addUniqueColumns(columns: {[tableName: string]: string}) {
    this.objHandlers.forEach((handler) => handler.addUniqueColumns(columns));
    return this;
  }

  public async addCustomSelectedColumns(columns: {[tableName: string]: string[]}) {
    this.objHandlers.forEach((handler) => handler.addCustomSelectedColumns(columns));
    return this;
  }

  public static async initFromConfig(config: DBExplorerConfig) {
    const exp = new DBExplorer(config.connectionName, config.schemaName);
    for (const [semType, entry] of Object.entries(config.entryPoints))
      await exp.addEntryPoint(semType, entry.table, entry.column);
    if (config.joinOptions)
      exp.addJoinOptions(config.joinOptions);
    if (config.headerNames)
      exp.addHeaderReplacers(config.headerNames);
    if (config.uniqueColumns)
      exp.addUniqueColumns(config.uniqueColumns);
    if (config.customSelectedColumns)
      await exp.addCustomSelectedColumns(config.customSelectedColumns);
    return exp;
  }

  public static async initFromConfigPath(_package: DG.Package, configPath: string = 'db-explorer/db-explorer-config.json') {
    try {
      const config = await _package.files.readAsText(configPath);
      return DBExplorer.initFromConfig(JSON.parse(config));
    } catch (e) {
      grok.shell.error(`Failed to load db-explorer config from ${configPath}`);
      console.error(e);
    }
    return null;
  }
}
