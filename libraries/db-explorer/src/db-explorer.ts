/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {DBExplorerConfig, EntityOptions, EntryPointOptions,
  QueryJoinOptions, ReferencedByObject, ReferenceObject, SchemaAndConnection} from './types';
import {DBExplorerObjectHandler, SemValueObjectHandler} from './object-handlers';
import {getDefaultRendererByName} from './renderer';
import {getSchemaReferences} from './utils';
import {identifierCode, TableNames} from './naming';

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
  /** What the Links section of every card of this explorer calls each table */
  public readonly names = new TableNames();
  constructor(
    private connectionName: string,
    private schemaName: string,
    private nqName?: string,
    private dataSourceName?: string,
    private catalogName?: string
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
      this.schemaName,
      this.catalogName,
      this.names,
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
    // a caller that waited for a load that failed gets no schema, as a later one does: with the connection alone it
    // would query a database that did not answer, and wait for that query to time out too
    return this._dbLoadPromise
      .then(() => this.loadingFailed ? null :
        ({schema: {references: this.references, referencedBy: this.referencedBy}, connection: this.connection} as SchemaAndConnection))
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

      // for implicit references, schema is always the primary one
      const schema = await getSchemaReferences(this.connection, this.schemaName, this.catalogName);
      this.references = schema.references;
      this.referencedBy = schema.referencedBy;
      this.schemasLoaded = true;
    } catch (e) {
      this.loadingFailed = true;
      // a database the stand cannot reach only means its object handlers are not offered
      console.warn('Failed to load DB schema, Object handlers not registered', e);
    } finally {
      if (this.connection)
        this.objHandlers.forEach((handler) => handler.connectionNqName = this.connection!.nqName); // set for detection in isApplicable
      // also when the connection is not found: the callers waiting for this load are answered
      this._dbLoadPromiseResolver();
    }
  }

  /** @param labels - what the Links section calls the row pointed at (label) and the rows pointing at it (reverseLabel) */
  public async addCustomRelation(tableName: string, columnName: string, refTable: string, refColumn: string, schemas?: {tableSchema?: string; refTableSchema?: string},
    labels?: {label?: string; reverseLabel?: string}) {
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
    schemaRefs[tableName][columnName] = {refTable, refColumn, refSchema: refSchemaName, label: labels?.label};
    if (!schemaRefBy[refTable])
      schemaRefBy[refTable] = {};
    schemaRefBy[refTable][refColumn] ??= [];
    // there already might be this reference, so first check if such ref already exists
    const existing = schemaRefBy[refTable][refColumn].find((r) => r.refTable === tableName && r.refColumn === columnName && r.refSchema === schemaName);
    if (existing) {
      existing.reverseLabel = labels?.reverseLabel ?? existing.reverseLabel;
      return;
    }
    schemaRefBy[refTable][refColumn].push({refTable: tableName, refColumn: columnName, refSchema: schemaName, reverseLabel: labels?.reverseLabel});
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
      this.nqName ?? this.connectionName,
      this.catalogName,
      this.names,
    );
    this.names.addEntryPoint({semType: semanticType, table: tableName, column: columnName,
      code: identifierCode(options?.regexpExample?.nonVariablePart, options?.matchRegexp)});
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

  /** Names and display of tables in the Links section of the cards, by table or schema.table */
  public addEntities(entities: {[tableName: string]: EntityOptions}) {
    this.names.addEntities(entities);
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
    const exp = new DBExplorer(config.connectionName, config.schemaName, config.nqName, config.dataSourceName, config.catalogName);
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
    if (config.entities)
      exp.addEntities(config.entities);
    if (config.explicitReferences) {
      config.explicitReferences.forEach((ref) => {
        exp.addCustomRelation(
          ref.table,
          ref.column,
          ref.refTable,
          ref.refColumn,
          {tableSchema: ref.schema, refTableSchema: ref.refSchema},
          {label: ref.label, reverseLabel: ref.reverseLabel},
        );
      });
    }
    if (config.customRenderers) {
      config.customRenderers.forEach((r) => {
        const rendererFunc = getDefaultRendererByName(r.renderer);
        const checkFunc = (tableName: string, colName: string, _value: string | number) => {
          return tableName?.toLowerCase() === r.table?.toLowerCase() && colName?.toLowerCase() === r.column?.toLowerCase();
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
