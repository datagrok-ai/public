/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {DBExplorerConfig, EntryPointOptions,
  QueryJoinOptions, ReferencedByObject, ReferenceObject, SchemaAndConnection} from './types';
import {DBExplorerObjectHandler, SemValueObjectHandler} from './object-handlers';
import {getDefaultRendererByName} from './renderer';

export class DBExplorer {
  private schemasLoaded: boolean = false;
  private connection: DG.DataConnection | null = null;
  private references: ReferenceObject = {};
  private referencedBy: ReferencedByObject = {};
  private _dbLoadPromise: Promise<void>;
  private _dbLoadPromiseResolver: (() => void) = () => {};
  private objHandlers: DBExplorerObjectHandler[] = [];
  private loadingFailed: boolean = false;
  public genericValueHandler: DBExplorerObjectHandler;
  constructor(
    private connectionName: string,
    private schemaName: string,
    private nqName?: string,
    private dataSourceName?: string
  ) {
    // in order for other functionality to be able to await the _dbLoadPromise,
    // we need to assign it to some promise and register resolver
    this._dbLoadPromise = new Promise<void>((resolve) => {
      this._dbLoadPromiseResolver = resolve;
    });
    this.genericValueHandler = new DBExplorerObjectHandler(
      {valueConverter: (a) => a,
        joinOptions: []},
      async () => {
        return await this.dbSchema;
      },
      nqName ?? connectionName,
      this.schemaName
    );
    DG.ObjectHandler.register(this.genericValueHandler);
    this.objHandlers.push(this.genericValueHandler);
  }

  private _dbLoadingPromise: Promise<void> | null = null;
  public get dbSchema(): Promise<SchemaAndConnection | null> {
    if (this.loadingFailed)
      return Promise.resolve(null);
    // to make sure it is only loaded once
    this._dbLoadingPromise ??= this.loadDbSchema();
    return this._dbLoadPromise
      .then(() => ({schema: {references: this.references, referencedBy: this.referencedBy}, connection: this.connection} as SchemaAndConnection))
      .catch((e) => {
        console.error(e);
        return null;
      });
  }

  /**
   * Note: this method is and should only be called once per instance
   * also, it is called automatically when dbSchema is requested for the first time
   * this method is lazy -- it starts loading only when someone requests it.
   */
  private async loadDbSchema() {
    try {
      const nqNameSplit = this.nqName?.split(':');
      const connections = nqNameSplit?.length === 2 ?
        await grok.dapi.connections.filter(`namespace = "${nqNameSplit[0]}:" and shortName = "${nqNameSplit[1]}"`).list() :
        await grok.dapi.connections.filter(`name="${this.connectionName}" or shortName = "${this.connectionName}"`).list();

      this.connection = connections.find((c) => (!this.nqName || c.nqName?.toLowerCase() === this.nqName.toLowerCase()) && (!this.dataSourceName || c.dataSource?.toLowerCase() === this.dataSourceName.toLowerCase())) ?? null;

      if (this.connection == null) {
        console.warn(`Connection ${this.connectionName} not found, Object handlers not registered`);
        this.loadingFailed = true;
        return;
      }

      const tables = await grok.dapi.connections.getSchema(this.connection, this.schemaName);
      if (!tables)
        throw new Error(`Schema ${this.schemaName} not found`);
      // for implicit references, schema is always the primary one
      this.references[this.schemaName] = {};
      this.referencedBy[this.schemaName] = {};
      const schemaRefs = this.references[this.schemaName];
      const schemaRefBy = this.referencedBy[this.schemaName];
      tables.forEach((table) => {
        const tableName = table.friendlyName ?? table.name;
        schemaRefs[tableName] = {};
        const t = schemaRefs[tableName];
        table.columns.forEach((column) => {
          const ref = column.referenceInfo;
          if (ref && ref.table && ref.column) {
            t[column.name] = {refTable: ref.table, refColumn: ref.column, refSchema: this.schemaName};
            if (!schemaRefBy[ref.table])
              schemaRefBy[ref.table] = {};

            if (!schemaRefBy[ref.table][ref.column])
              schemaRefBy[ref.table][ref.column] = [];

            schemaRefBy[ref.table][ref.column].push({refTable: tableName, refColumn: column.name, refSchema: this.schemaName});
          }
        });
      });
      this.schemasLoaded = true;
    } catch (_e) {
      this.loadingFailed = true;
      console.warn('Failed to load DB schema, Object handlers not registered');
      console.error(_e);
    }
    if (this.connection)
      this.objHandlers.forEach((handler) => handler.connectionNqName = this.connection!.nqName); // set for detection in isApplicable
    this._dbLoadPromiseResolver();
  }

  public async addCustomRelation(tableName: string, columnName: string, refTable: string, refColumn: string, schemas?: {tableSchema?: string; refTableSchema?: string}) {
    if (!this.schemasLoaded) {
      try {
        await this._dbLoadPromise;
      } catch (_e) {
        console.error(_e);
      }
    }
    if (!this.schemasLoaded || !this.connection || !this.references || !this.referencedBy) {
      console.warn('Failed to add custom relation, DB schema not loaded');
      return;
    }
    const schemaName = schemas?.tableSchema ?? this.schemaName;
    const refSchemaName = schemas?.refTableSchema ?? this.schemaName;
    if (!this.references[schemaName])
      this.references[schemaName] = {};
    if (!this.referencedBy[refSchemaName])
      this.referencedBy[refSchemaName] = {};
    const schemaRefs = this.references[schemaName];
    const schemaRefBy = this.referencedBy[refSchemaName];
    if (!schemaRefs[tableName])
      schemaRefs[tableName] = {};
    schemaRefs[tableName][columnName] = {refTable, refColumn, refSchema: refSchemaName};
    if (!schemaRefBy[refTable])
      schemaRefBy[refTable] = {};
    schemaRefBy[refTable][refColumn] ??= [];
    schemaRefBy[refTable][refColumn].push({refTable: tableName, refColumn: columnName, refSchema: schemaName});
  }

  public addEntryPoint(
    semanticType: string,
    tableName: string,
    columnName: string,
    options?: Partial<EntryPointOptions>
  ) {
    const schemaPromise = () => this.dbSchema;
    const fullOpts = {...{valueConverter: (a: string | number) => a, joinOptions: []}, ...(options ?? {})};
    const handler = new SemValueObjectHandler(
      semanticType,
      tableName,
      columnName,
      fullOpts,
      schemaPromise,
      this.schemaName,
      this.nqName ?? this.connectionName
    );
    // if the matching regexp is also provided, register it
    if (options?.matchRegexp)
      DG.SemanticValue.registerRegExpDetector(semanticType, options.matchRegexp);
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

  public addCustomSelectedColumns(columns: {[tableName: string]: string[]}) {
    this.objHandlers.forEach((handler) => handler.addCustomSelectedColumns(columns));
    return this;
  }

  public static initFromConfig(config: DBExplorerConfig) {
    const exp = new DBExplorer(config.connectionName, config.schemaName, config.nqName, config.dataSourceName);
    for (const [semType, entry] of Object.entries(config.entryPoints))
      exp.addEntryPoint(semType, entry.table, entry.column, {regexpExample: entry.regexpExample, matchRegexp: entry.matchRegexp});
    if (config.joinOptions)
      exp.addJoinOptions(config.joinOptions);
    if (config.headerNames)
      exp.addHeaderReplacers(config.headerNames);
    if (config.uniqueColumns)
      exp.addUniqueColumns(config.uniqueColumns);
    if (config.customSelectedColumns)
      exp.addCustomSelectedColumns(config.customSelectedColumns);
    if (config.explicitReferences) {
      config.explicitReferences.forEach((ref) => {
        exp.addCustomRelation(
          ref.table,
          ref.column,
          ref.refTable,
          ref.refColumn,
          {tableSchema: ref.schema, refTableSchema: ref.refSchema}
        );
      });
    }
    if (config.customRenderers) {
      config.customRenderers.forEach((r) => {
        const rendererFunc = getDefaultRendererByName(r.renderer);
        const checkFunc = (tableName: string, colName: string, _value: string | number) => {
          return tableName === r.table && colName === r.column;
        };
        exp.addCustomRenderer(checkFunc, rendererFunc);
      });
    }
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
