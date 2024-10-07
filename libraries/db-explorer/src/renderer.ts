/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {DBValueObject, EntryPointOptions, QueryJoinOptions, SchemaInfo} from './types';
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

export function imageRenderer(fullUrl: string) {
  const nameHost = textRenderer(fullUrl);
  const loaderDiv = getLoaderDiv();
  const host = ui.divH([nameHost, loaderDiv]);
  ui.tooltip.bind(nameHost, 'Unable to get image');

  async function replaceWithImage() {
    try {
      const qRes = await grok.dapi.fetchProxy(fullUrl, {});
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

export function ownIdRenderer(id: string | number, tableName: string, colName: string) {
  const nameHost = textRenderer(id.toString(), false);
  nameHost.style.color = 'var(--blue-1)';
  nameHost.style.cursor = 'pointer';
  ui.tooltip.bind(nameHost, 'Click to make it current object');
  nameHost.addEventListener('click', (e) => {
    e.stopImmediatePropagation();
    grok.shell.o = new DBValueObject(tableName, colName, id);
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

  private headerReplacers: {[tableName: string]: string} = {};
  private uniqueColNames: {[tableName: string]: string} = {};
  private customSelectedColumns: {[tableName: string]: Set<string>} = {};
  private defaultHeaderReplacerColumns: string[] = ['name'];
  constructor(
    private schemaInfo: SchemaInfo,
    private connection: DG.DataConnection,
    private schemaName: string,
    private entryPointOptions: EntryPointOptions
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

          ui.tooltip.bind(nameHost, () => this.renderDataFrame(clearedDF, tableName));
          nameHost.addEventListener('click', (e) => {
            e.stopImmediatePropagation();
            grok.shell.o = new DBValueObject(tableName, colName, id);
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
    const res = await queryDB(this.connection, tableName, match, matchValue, this.schemaName, joinOptions);
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
        return this.renderDataFrame(rowDf, tableName);
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
    return this.renderDataFrame(df, tableName);
  }

  renderDataFrame(df: DG.DataFrame, tableName: string) {
    const clearedDF = removeEmptyCols(df);
    let entries = Object.entries(clearedDF.toJson()[0]) as [string, string | number][];

    if (this.customSelectedColumns[tableName])
      entries = entries.filter(([key]) => this.customSelectedColumns[tableName].has(key));

    const mainTable = ui.table(entries, (entry) => {
      const key = entry[0];
      const value = entry[1];
      const cutomRenderer = this.customRenderers.find((renderer) => renderer.check(tableName, key, value));
      if (cutomRenderer) return [textRenderer(key), cutomRenderer.renderer(value, this.connection)];

      const refInfo = this.schemaInfo.references?.[tableName]?.[key];
      if (refInfo)
        return [textRenderer(key), this.refIdRenderer(this.connection, value, refInfo.refTable, refInfo.refColumn)];
      const isUnqueCol = this.uniqueColNames[tableName] === key;
      if (isUnqueCol)
        return [textRenderer(key), ownIdRenderer(value, tableName, key)];


      return [textRenderer(key), textRenderer(value)];
    });
    return ui.divV([mainTable]);
  }

  renderAssociations(acc: DG.Accordion, schema: SchemaInfo, curTable: string, curDf: DG.DataFrame) {
    const attachOpenInWorkspaceMenu = (tableName: string, colName: string, value: any, pane: DG.AccordionPane) => {
      const menu = DG.Menu.popup();
      menu.item(`Add all associated ${tableName} entries to workspace`, async () => {
        const pi = DG.TaskBarProgressIndicator.create('Opening all associated entries');
        const assocDf = await this.getTable(tableName, colName, value, this.entryPointOptions.joinOptions);
        if (assocDf) {
          assocDf.name = `${tableName}`;
          grok.shell.addTableView(assocDf);
        }
        pi.close();
      });
      pane.root.addEventListener('contextmenu', (e) => {
        e.preventDefault();
        setTimeout(() => menu.show());
      });
    };

    const refTable = schema.referencedBy[curTable];
    if (!refTable || Object.keys(refTable).length == 0) return;
    const _linksPane = acc.addPane('Links', () => {
      const linksAcc = ui.accordion(`Links to ${curTable}`);
      Object.entries(refTable).forEach(([refedColumn, refInfo]) => {
        const val = curDf.col(refedColumn)?.get(0);
        if (!val) return;
        const _pane = linksAcc.addPane(refedColumn, () => {
          const colAcc = ui.accordion(`Links to ${curTable}_${refedColumn}`);
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
                  if (res.rowCount > MAX_MULTIROW_VALUES)
                    singleColPane.name = `${singleColPane.name} (${MAX_MULTIROW_VALUES} / ${res.rowCount})`;
                  return res.root;
                });
              }
              const multiTableAcc = ui.accordion(`Multiple links to ${refTable}`);
              refInfos.forEach((ref) => {
                const colPane = multiTableAcc.addPane(ref.refColumn, () => {
                  return ui.wait(async () => {
                    const res = await this.renderMultiRowTable(ref.refTable, ref.refColumn, val, this.entryPointOptions.joinOptions);
                    if (res.rowCount > MAX_MULTIROW_VALUES)
                      colPane.name = `${singleColPane.name} (${MAX_MULTIROW_VALUES} / ${res.rowCount})`;
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

      return linksAcc.root;
    });
  }
}
