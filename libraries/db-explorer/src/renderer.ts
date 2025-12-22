/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {DBExplorerConfig, DBValueObject, EntryPointOptions, QueryJoinOptions, ReferencedByObject, SchemaAndConnection, SupportedRenderer} from './types';
import {queryDB} from './query';

export const MAX_MULTIROW_VALUES = 10;
/** normal renderer supporting copy and elipsis */
export function textRenderer(value: string | number, withTooltip = true) {
  const nameHost = ui.div(value.toString());
  nameHost.style.maxWidth = '200px';
  nameHost.style.overflow = 'hidden';
  nameHost.style.whiteSpace = 'nowrap';
  nameHost.style.textOverflow = 'ellipsis';

  // ability to copy
  const menu = DG.Menu.popup();
  menu.item('Copy', () => {
    navigator?.clipboard?.writeText(value.toString());
  });
  nameHost.addEventListener('contextmenu', (e) => {
    e.preventDefault();
    e.stopImmediatePropagation();
    setTimeout(() => menu.show());
  });

  // approximately what will fit in 150 px
  if (value.toString().length > 20 && withTooltip) ui.tooltip.bind(nameHost, value.toString());
  return nameHost;
}

export function moleculeRenderer(value: string) {
  try {
    const molDiv = grok.chem.drawMolecule(value, 200, 200);
    return molDiv;
  } catch (e) {
    console.error(e);
  }

  return textRenderer(value);
}

export function imageRenderer(fullUrl: string, useProxy: boolean = true) {
  const nameHost = textRenderer(fullUrl);
  const loaderDiv = getLoaderDiv();
  const host = ui.divH([nameHost, loaderDiv]);
  ui.tooltip.bind(nameHost, 'Unable to get image');

  async function replaceWithImage() {
    try {
      const qRes = useProxy ? await grok.dapi.fetchProxy(fullUrl, {}) : await fetch(fullUrl);
      if (!qRes) throw new Error('');
      const blob = await qRes.blob();
      if (!blob) throw new Error('');
      const canvas = ui.canvas(200, 200);
      const ctx = canvas.getContext('2d');
      if (!ctx) throw new Error('');
      const image = new Image();
      image.onload = () => {
        ctx.drawImage(image, 0, 0, 200, 200);
      };
      image.src = URL.createObjectURL(blob);

      nameHost.remove();
      host.appendChild(canvas);
    } catch (e) {
      console.error(e);
    } finally {
      loaderDiv.remove();
    }
  }
  replaceWithImage();
  return host;
}

export function helmRenderer(value: string): HTMLElement {
  return ui.wait(async () => {
    //@ts-ignore
    const helmInput = await ui.input.helmAsync('helm', {
      editable: false,
    });
    helmInput.setStringValue(value);
    await DG.delay(200); // wait for proper sizing
    helmInput.getInput().addEventListener('click', () => {
      grok.shell.o = helmInput.getValue();
    });
    helmInput.getInput().addEventListener('dblclick', () => {
      helmInput.showEditorDialog();
    });

    helmInput.getInput().style.width = '100%';
    helmInput.getInput().style.setProperty('height', '300px', 'important');
    return helmInput.getInput();
  });
}

export function rawImageRenderer(rawImage: string) {
  const nameHost = textRenderer('image');
  const loaderDiv = getLoaderDiv();
  const host = ui.divH([nameHost, loaderDiv]);
  ui.tooltip.bind(nameHost, 'Unable to render image');
  async function replaceWithImage() {
    try {
      if (!rawImage) throw new Error('Empty image string');
      const canvas = ui.canvas(200, 200);
      const ctx = canvas.getContext('2d');
      if (!ctx) throw new Error('');
      const image = new Image();
      image.onload = () => {
        ctx.drawImage(image, 0, 0, 200, 200);
      };
      image.src = 'data:image/png;base64,' + rawImage;

      nameHost.remove();
      host.appendChild(canvas);
    } catch (e) {
      console.error(e);
    } finally {
      loaderDiv.remove();
    }
  }
  replaceWithImage();
  return host;
}

export function ownIdRenderer(id: string | number, tableName: string, colName: string, connectionNqName: string, schemaName: string) {
  const nameHost = textRenderer(id.toString(), false);
  nameHost.style.color = 'var(--blue-1)';
  nameHost.style.cursor = 'pointer';
  ui.tooltip.bind(nameHost, 'Click to explore');
  nameHost.addEventListener('click', (e) => {
    e.stopImmediatePropagation();
    grok.shell.o = new DBValueObject(connectionNqName, schemaName, tableName, colName, id);
  });
  return nameHost;
}

function removeEmptyCols(df: DG.DataFrame) {
  df.columns.names().forEach((colName) => {
    if (!df.col(colName) || df.col(colName)!.isNone(0)) df.columns.remove(colName);
  });
  return df;
}

export function getLoaderDiv() {
  const div = ui.div([], {style: {width: '50px', height: '24px', position: 'relative'}});
  div.innerHTML = `<div class="grok-loader"><div></div><div></div><div></div><div></div></div>`;
  return div;
}

export class DBExplorerRenderer {
  private valueReplacers: {
    check: (tableName: string, colName: string, value: string | number) => boolean;
    replacer: (df: DG.DataFrame, value: string | number) => string;
  }[] = [];
  private customRenderers: {
    check: (tableName: string, colName: string, value: string | number) => boolean;
    renderer: (value: string | number, connection: DG.DataConnection) => HTMLElement;
  }[] = [];

  private semtypeRenderers = [{
    check: (col?: DG.Column) => col?.semType === DG.SEMTYPE.MOLECULE,
    renderer: (value: any) => moleculeRenderer(value),
  }, {
    check: (col?: DG.Column) => col?.semType === DG.SEMTYPE.MACROMOLECULE && col?.meta?.units === 'helm',
    renderer: (value: any) => helmRenderer(value),
  },
  {
    check: (col?: DG.Column) => col?.semType === 'rawPng',
    renderer: (value: any) => rawImageRenderer(value),
  },
  {
    check: (col?: DG.Column) => col?.semType === 'ImageUrl',
    renderer: (value: any) => imageRenderer(value),
  }
  ];

  private headerReplacers: {[tableName: string]: string} = {};
  private uniqueColNames: {[tableName: string]: string} = {};
  private customSelectedColumns: {[tableName: string]: Set<string>} = {};
  private defaultHeaderReplacerColumns: string[] = ['name', 'type', 'standard_type'];
  constructor(
    protected schemaInfoPromise: () => Promise<SchemaAndConnection | null>,
    protected schemaName: string,
    protected entryPointOptions: EntryPointOptions
  ) {}

  addHeaderReplacers(replacers: {[tableName: string]: string}) {
    this.headerReplacers = {...this.headerReplacers, ...replacers};
  }

  addDefaultHeaderReplacerColumns(columns: string[]) {
    this.defaultHeaderReplacerColumns = [...this.defaultHeaderReplacerColumns, ...columns];
  }

  addCustomSelectedColumns(columns: {[tableName: string]: string[]}) {
    const columnSets: {[tableName: string]: Set<string>} = {};
    Object.entries(columns).forEach(([tableName, cols]) => {
      columnSets[tableName] = new Set(cols);
    });
    this.customSelectedColumns = {...this.customSelectedColumns, ...columnSets};
  }

  addUniqueColNames(cols: {[tableName: string]: string}) {
    this.uniqueColNames = {...this.uniqueColNames, ...cols};
  }

  addValueReplacer(
    check: (tableName: string, colName: string, value: string | number) => boolean,
    replacer: (df: DG.DataFrame, value: string | number) => string
  ) {
    this.valueReplacers.push({check, replacer});
  }
  addCustomRenderer(
    check: (tableName: string, colName: string, value: string | number) => boolean,
    renderer: (value: string | number, connection: DG.DataConnection) => HTMLElement
  ) {
    this.customRenderers.push({check, renderer});
  }

  refIdRenderer(connection: DG.DataConnection, id: string | number, tableName: string, colName: string) {
    const loaderDiv = getLoaderDiv();
    const nameHost = textRenderer(id.toString(), false);
    const host = ui.divH([nameHost, loaderDiv]);
    function unableToRetrieve(reason: string = `Unable to retrieve ${tableName} information`) {
      loaderDiv.remove();
      ui.tooltip.bind(nameHost, reason);
    }
    if (id == null || id === '') {
      unableToRetrieve('ID is empty');
      return host;
    }

    try {
      queryDB(connection, tableName, colName, id, this.schemaName, this.entryPointOptions.joinOptions)
        .then((df) => {
          if (df.rowCount == 0) {
            unableToRetrieve();
            return;
          }
          const clearedDF = removeEmptyCols(df);
          const replacer = this.valueReplacers.find((replacer) => replacer.check(tableName, colName, id));
          if (replacer) nameHost.textContent = replacer.replacer(clearedDF, id);

          ui.tooltip.bind(nameHost, () => ui.wait(async () => await this.renderDataFrame(clearedDF, tableName)));
          nameHost.addEventListener('click', (e) => {
            e.stopImmediatePropagation();
            this.schemaInfoPromise().then((schemaAndConnection) => {
              if (schemaAndConnection && schemaAndConnection.connection && schemaAndConnection.schema)
                grok.shell.o = new DBValueObject(schemaAndConnection.connection.nqName, this.schemaName, tableName, colName, id);
            });
          });
          nameHost.style.color = 'var(--blue-1)';
          nameHost.style.cursor = 'pointer';
        })
        .finally(() => {
          loaderDiv.remove();
        });
    } catch (e) {
      console.error(e);
      unableToRetrieve();
    }
    return host;
  }

  async getTable(tableName: string, match: string, matchValue: string | number, joinOptions: QueryJoinOptions[] = []) {
    const schemaAndConnection = await this.schemaInfoPromise();
    const res = await queryDB(schemaAndConnection?.connection ?? null, tableName, match, matchValue, this.schemaName, joinOptions);
    return res;
  }

  async renderMultiRowTable(
    tableName: string,
    match: string,
    matchValue: string | number,
    joinOptions: QueryJoinOptions[] = []
  ) {
    const df = await this.getTable(tableName, match, matchValue, joinOptions);
    if (df.rowCount < 1)
      return {root: ui.divText('ID not found'), rowCount: 0};
    if (df.rowCount === 1)
      return {root: this.renderDataFrame(df, tableName), rowCount: 1};
    const acc = ui.accordion(`Multiple rows for ${tableName}`);
    const dfMaxCount = Math.min(df.rowCount, 10);
    const replaceColName = this.headerReplacers[tableName];
    for (let i = 0; i < dfMaxCount; i++) {
      let paneName = `Row ${i + 1}`;
      if (replaceColName && df.col(replaceColName)?.get(i)) {
        paneName = df.col(replaceColName)!.get(i).toString();
      } else {
        const f = this.defaultHeaderReplacerColumns.find((col) => df.col(col) && df.col(col)!.get(i));
        if (f) paneName = df.col(f)!.get(i).toString();
      }
      acc.addPane(paneName, () => {
        const rowBitset = DG.BitSet.create(df.rowCount);
        rowBitset.set(i, true);
        const rowDf = df.clone(rowBitset);
        return ui.wait(async () => await this.renderDataFrame(rowDf, tableName));
      });
    }
    return {root: acc.root, rowCount: df.rowCount};
  }

  async renderTable(
    tableName: string,
    match: string,
    matchValue: string | number,
    joinOptions: QueryJoinOptions[] = []
  ) {
    const df = await this.getTable(tableName, match, matchValue, joinOptions);
    if (df.rowCount < 1) return ui.divText('ID not found');
    await df.meta.detectSemanticTypes();
    return this.renderDataFrame(df, tableName);
  }

  /**
   * Renders dataframe as a card with all associated references clickable
   * In cases where your table is aggregated from different db tables, you can pass columnTablesMap property to specify
   * which column belongs to which table in the database schema
   * @param df
   * @param tableName
   * @param options
   * @returns
   */
  async renderDataFrame(df: DG.DataFrame, tableName: string, options?: {keepEmptyValues: boolean, skipCustomSelected?: boolean}) {
    const schemaAndConnection = await this.schemaInfoPromise();
    if (!schemaAndConnection)
      return ui.divText('Schema information is not available');
    if (df.rowCount === 0)
      return ui.divText('No data');

    const clearedDF = options?.keepEmptyValues ? df: removeEmptyCols(df);
    let entries = clearedDF.columns.names()
      .map((colName) => [colName, clearedDF.col(colName)!.isNone(0) ? '' : clearedDF.col(colName)!.getString(0)]) as [string, string | number][];

    if (!options?.skipCustomSelected && this.customSelectedColumns[tableName]) {
      entries = entries.filter(([key]) => this.customSelectedColumns[tableName].has(key));
      // reorder entries according to customSelectedColumns if present
      const cols = Array.from(this.customSelectedColumns[tableName]);
      entries.sort((a, b) => {
        return cols.indexOf(a[0]) - cols.indexOf(b[0]);
      });
    }
    // similarly here, the column name in the df might not match the db column name if the df is a result of a query joining multiple tables
    const dbColumnNames: {[colName: string]: string} = {};
    const dbColTableNames: {[colName: string]: string} = {};
    entries.forEach(([colName, _]) => {
      const col = clearedDF.col(colName);
      if (col?.getTag(DG.Tags.DbColumn))
        dbColumnNames[colName] = col.getTag(DG.Tags.DbColumn)!;
      if (col?.getTag(DG.Tags.DbTable))
        dbColTableNames[colName] = col.getTag(DG.Tags.DbTable)!;
    });

    const isAggregatedFromMultipleTables = Object.keys(dbColTableNames).some((colName) => dbColTableNames[colName] !== tableName);

    const mainTable = ui.table(entries, (entry) => {
      const colName = entry[0];
      const dbColName = dbColumnNames[colName] ?? colName;
      const value = entry[1];
      const dbTableName = dbColTableNames[colName] ?? tableName; // here, it is mapped from the df column, not its db name
      const cutomRenderer = this.customRenderers.find((renderer) => renderer.check(dbTableName, dbColName, value));
      if (cutomRenderer) return [textRenderer(colName), cutomRenderer.renderer(value, schemaAndConnection.connection)];
      // if its autodetected as a molecule or helm or something else in the future
      const semTypeRenderer = this.semtypeRenderers.find((renderer) => {
        const col = clearedDF.col(colName);
        return renderer.check(col ?? undefined);
      });
      if (semTypeRenderer) return [textRenderer(colName), semTypeRenderer.renderer(value)];

      const refInfo = schemaAndConnection.schema.references?.[dbTableName]?.[dbColName];
      if (refInfo)
        return [textRenderer(colName), this.refIdRenderer(schemaAndConnection.connection, value, refInfo.refTable, refInfo.refColumn)];
      // if the given column is referenced by other tables, then also use refIdRenderer
      const refedByInfo = schemaAndConnection.schema.referencedBy?.[dbTableName]?.[dbColName];
      if (refedByInfo && refedByInfo.length > 0)
        return [textRenderer(colName), this.refIdRenderer(schemaAndConnection.connection, value, dbTableName, dbColName)];

      const isUnqueCol = this.uniqueColNames[dbTableName] === dbColName;
      if (isUnqueCol) {
        if (!isAggregatedFromMultipleTables) // if its not aggregated, the tooltip in of the id will be the exact same as card, so we use simpler ownerIdRenderer
          return [textRenderer(colName), ownIdRenderer(value, dbTableName, dbColName, schemaAndConnection.connection.nqName, this.schemaName)];
        return [textRenderer(colName), this.refIdRenderer(schemaAndConnection.connection, value, dbTableName, dbColName)]; // otherwise use refIdRenderer
      }


      return [textRenderer(colName), textRenderer(value)];
    });
    return ui.divV([mainTable]);
  }

  async renderAssociations(acc: DG.Accordion, schemaPromise: () => Promise<SchemaAndConnection | null>, curTable: string, curDf: DG.DataFrame) {
    const schemaAndConnection = await schemaPromise();
    if (!schemaAndConnection)
      return;
    const addAllAssociated = async (tableName: string, colName: string, value: any) => {
      const pi = DG.TaskBarProgressIndicator.create('Opening all associated entries');
      try {
        const assocDf = await this.getTable(tableName, colName, value, this.entryPointOptions.joinOptions);
        if (assocDf) {
          if (assocDf.rowCount === 0) {
            grok.shell.info(`No associated ${tableName} entries found`);
            return;
          }
          assocDf.name = `${tableName}`;
          grok.shell.addTableView(assocDf);
        }
      } catch (e) {
        console.error(e);
        grok.shell.error(`Failed to open associated ${tableName} entries`);
      } finally {
        pi.close();
      }
    };

    const addIconToPane = (pane: DG.AccordionPane, tableName: string, colName: string, value: any) => {
      const icon = ui.icons.add(() => {}, `Add all associated ${tableName} entries to workspace`);
      // need separate event to avoid click on accordion pane
      icon.addEventListener('click', (e) => {
        e.stopImmediatePropagation();
        e.preventDefault();
        addAllAssociated(tableName, colName, value);
      });
      pane.root.getElementsByClassName('d4-accordion-pane-header')?.[0]?.appendChild(icon);
    };

    const attachOpenInWorkspaceMenu = (tableName: string, colName: string, value: any, pane: DG.AccordionPane) => {
      const menu = DG.Menu.popup();


      menu.item(`Add all associated ${tableName} entries to workspace`, async () => {
        addAllAssociated(tableName, colName, value);
      });
      pane.root.addEventListener('contextmenu', (e) => {
        e.preventDefault();
        setTimeout(() => menu.show());
      });
      addIconToPane(pane, tableName, colName, value);
    };
    // columns in the passed dataframe might be from different tables in the db schema (as a result of queries)
    // depending on from where this dataFrame come from

    const colList = curDf.columns.names();
    const dbTableNameSet = new Set<string>();
    dbTableNameSet.add(curTable);
    // build map of actual column name to db column name
    const dbNameToActColNameMap: {[tableName: string]: {[colName: string]: string}} = {};
    colList.forEach((colName) => {
      const dbTableName = curDf.col(colName)!.getTag(DG.Tags.DbTable);
      if (dbTableName)
        dbTableNameSet.add(dbTableName);
      const dbColName = curDf.col(colName)!.getTag(DG.Tags.DbColumn);
      if (dbColName) {
        const actTableName = dbTableName ?? curTable;
        if (!dbNameToActColNameMap[actTableName])
          dbNameToActColNameMap[actTableName] = {};
        dbNameToActColNameMap[actTableName][dbColName] = colName;
      }
    });
    // all references to all tables referencing any of the db tables in the current dataframe
    const refTables = Array.from(dbTableNameSet).map((t) => [t, schemaAndConnection.schema.referencedBy[t]] as [string, ReferencedByObject[string]]).filter((t) => t[1] != null && Object.keys(t[1]).length > 0);
    if (refTables.length === 0)
      return;

    const _linksPane = acc.addPane('Links', () => {
      const linksAcc = ui.accordion(`Links to ${curTable}`);
      if ((curDf?.rowCount ?? 0) === 0)
        return ui.divText('No data');
      refTables.forEach(([actTableName, refTable]) => {
        const prefix = refTables.length > 1 ? `${actTableName}.` : '';
        Object.entries(refTable).forEach(([refedColumn, refInfo]) => {
          const tableColActName = dbNameToActColNameMap[actTableName][refedColumn] ?? refedColumn;
          const val = curDf.col(tableColActName)?.get(0);
          if (!val) return;
          const _pane = linksAcc.addPane(prefix + refedColumn, () => {
            const colAcc = ui.accordion(`Links to ${actTableName}_${refedColumn}`);
            // there can be cases when same column is referneced by two or more columns in other table.
            // example is chembl table molecule_hiearchy where all 3 columns are referencing to molregno
            // need to account for such cases
            const tableRefs = new Map<string, { refTable: string; refColumn: string }[]>();
            refInfo.forEach((ref) => {
              if (!tableRefs.has(ref.refTable)) tableRefs.set(ref.refTable, []);
            tableRefs.get(ref.refTable)!.push(ref);
            });

            tableRefs.forEach((refInfos, refTable) => {
              const singleColPane = colAcc.addPane(refTable, () => {
                if (refInfos.length === 1) {
                  return ui.wait(async () => {
                    const res = await this.renderMultiRowTable(refInfos[0].refTable, refInfos[0].refColumn, val, this.entryPointOptions.joinOptions);
                    if (res.rowCount > MAX_MULTIROW_VALUES) {
                      singleColPane.name = `${singleColPane.name} (${MAX_MULTIROW_VALUES} / ${res.rowCount})`;
                      addIconToPane(singleColPane, refInfos[0].refTable, refInfos[0].refColumn, val);
                    }
                    return res.root;
                  });
                }
                const multiTableAcc = ui.accordion(`Multiple links to ${refTable}`);
                refInfos.forEach((ref) => {
                  const colPane = multiTableAcc.addPane(ref.refColumn, () => {
                    return ui.wait(async () => {
                      const res = await this.renderMultiRowTable(ref.refTable, ref.refColumn, val, this.entryPointOptions.joinOptions);
                      if (res.rowCount > MAX_MULTIROW_VALUES) {
                        colPane.name = `${singleColPane.name} (${MAX_MULTIROW_VALUES} / ${res.rowCount})`;
                        addIconToPane(colPane, ref.refTable, ref.refColumn, val);
                      }
                      return res.root;
                    });
                  });
                  attachOpenInWorkspaceMenu(ref.refTable, ref.refColumn, val, colPane);
                });
                return multiTableAcc.root;
              });
              if (refInfos.length === 1)
                attachOpenInWorkspaceMenu(refInfos[0].refTable, refInfos[0].refColumn, val, singleColPane);
            });
            return colAcc.root;
          });
        });
      });

      return linksAcc.root;
    });
  }
}


export class ExampleExplorerRenderer extends DBExplorerRenderer {
  private static schemaCache: Map<string, Promise<SchemaAndConnection | null>> = new Map();

  constructor(
    connection: DG.DataConnection,
    shemaName: string,
    entryPointOptions: EntryPointOptions
  ) {
    const cacheKey = `${connection.id}||${shemaName}`;
    if (ExampleExplorerRenderer.schemaCache.has(cacheKey)) {
      const schemaPromise = ExampleExplorerRenderer.schemaCache.get(cacheKey)!;
      super(() => schemaPromise, shemaName, entryPointOptions);
    } else {
      const schemaPromise: Promise<SchemaAndConnection> = grok.dapi.connections.getSchema(connection, shemaName).then((tables) => {
        const references: {[tableName: string]: {[colName: string]: {refTable: string; refColumn: string}}} = {};
        const referencedBy: {[tableName: string]: {[colName: string]: {refTable: string; refColumn: string}[]}} = {};
        tables.forEach((table) => {
          const tableName = table.friendlyName ?? table.name;

          references[tableName] = {};
          const t = references[tableName];
          table.columns.forEach((column) => {
            const ref = column.referenceInfo;
            if (ref && ref.table && ref.column) {
              t[column.name] = {refTable: ref.table, refColumn: ref.column};
              if (!referencedBy[ref.table])
                referencedBy[ref.table] = {};

              if (!referencedBy[ref.table][ref.column])
                referencedBy[ref.table][ref.column] = [];

              referencedBy[ref.table][ref.column].push({refTable: tableName, refColumn: column.name});
            }
          });
        });
        return {
          schema: {references, referencedBy},
          connection,
        };
      });
      ExampleExplorerRenderer.schemaCache.set(cacheKey, schemaPromise);
      super(() => schemaPromise, shemaName, entryPointOptions);
    }
  }

  override async renderDataFrame(df: DG.DataFrame, tableName: string, _opts?: {keepEmptyValues: boolean, columnTablesMap?: {[colName: string]: string}, skipCustomSelected?: boolean}): Promise<HTMLDivElement> {
    return super.renderDataFrame(df, tableName, {keepEmptyValues: true});
  }

  override async renderTable(tableName: string, match: string, matchValue: string | number, joinOptions?: QueryJoinOptions[]): Promise<HTMLDivElement> {
    const schemaAndConnection = await this.schemaInfoPromise();
    if (!schemaAndConnection)
      return ui.divText('Schema information is not available');
    const df = await queryDB(schemaAndConnection.connection, tableName, '', '', this.schemaName, joinOptions, true);
    if (df.rowCount < 1) return ui.divText('ID not found');
    await df.meta.detectSemanticTypes();
    return this.renderDataFrame(df, tableName);
  }
}

export function renderExampleCard(connection: DG.DataConnection, schemaName: string, tableName: string,
  config: Partial<DBExplorerConfig>): HTMLElement {
  const ex = new ExampleExplorerRenderer(
    connection,
    schemaName,
    {joinOptions: config.joinOptions ?? [], valueConverter: (a) => a}
  );
  ex.addCustomSelectedColumns(config.customSelectedColumns ?? {});
  (config.customRenderers ?? []).forEach((cr) =>
    ex.addCustomRenderer((t, c, _v) => t === cr.table && c === cr.column, (v) => getDefaultRendererByName(cr.renderer)(v))
  );
  (config.uniqueColumns) && ex.addUniqueColNames(config.uniqueColumns);
  (config.customSelectedColumns) && ex.addCustomSelectedColumns(config.customSelectedColumns);
  const c = ui.card(ui.wait(async () => {
    return await ex.renderTable(tableName, '', '', config.joinOptions ?? []);
  }));
  c.style.width = 'unset';
  return c;
}

export function getDefaultRendererByName(rendererName: SupportedRenderer): (value: string | number) => HTMLElement {
  switch (rendererName) {
    case 'rawImage':
      return (v) => rawImageRenderer(v as string);
    case 'imageURL':
      return (v) => imageRenderer(v as string, true);
    case 'molecule':
      return (v) => moleculeRenderer(v as string);
    case 'helm':
      return (v) => helmRenderer(v as string);
    default:
      return (v) => textRenderer(v as string, true);
  }
}
