/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {DBExplorerConfig, DBValueObject, EntryPointOptions, LinkGroup, LinkRef, QueryJoinOptions, SchemaAndConnection, SupportedRenderer} from './types';
import {queryDB} from './query';
import {getSchemaReferences} from './utils';
import {ciEquals, ciGet, lcFirst, linkGroups, roleWords, TableNames, withoutParent} from './naming';

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
    // the platform puts a JS input's element into the input's own root one timer tick after the input is made
    // (JsInputProxy, Timer.run): the drawing is taken out of it after that, or that move takes it back out of the card
    await new Promise((resolve) => setTimeout(resolve));
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

export function ownIdRenderer(id: string | number, tableName: string, colName: string, connectionNqName: string, schemaName: string, catalogName?: string) {
  const nameHost = textRenderer(id.toString(), false);
  nameHost.style.color = 'var(--blue-1)';
  nameHost.style.cursor = 'pointer';
  ui.tooltip.bind(nameHost, 'Click to explore');
  nameHost.addEventListener('click', (e) => {
    e.stopImmediatePropagation();
    grok.shell.o = new DBValueObject(connectionNqName, schemaName, tableName, colName, id, undefined, catalogName);
  });
  return nameHost;
}

function removeEmptyCols(df: DG.DataFrame) {
  const dfClone = df.clone();
  df.columns.names().forEach((colName) => {
    if (!df.col(colName) || df.col(colName)!.isNone(0)) dfClone.columns.remove(colName);
  });
  return dfClone;
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
    protected entryPointOptions: EntryPointOptions,
    protected catalogName?: string,
    protected names: TableNames = new TableNames(),
  ) {}

  addHeaderReplacers(replacers: {[tableName: string]: string}) {
    this.headerReplacers = {...this.headerReplacers, ...replacers};
  }

  addDefaultHeaderReplacerColumns(columns: string[]) {
    this.defaultHeaderReplacerColumns = [...this.defaultHeaderReplacerColumns, ...columns];
  }

  /** schema qualified table names and corresponding column names */
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

  refIdRenderer(connection: DG.DataConnection, id: string | number, schemaName: string, tableName: string, colName: string) {
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
      queryDB(connection, tableName, colName, id, schemaName, this.entryPointOptions.joinOptions, false, this.catalogName)
        .then((df) => {
          if (df.rowCount == 0) {
            unableToRetrieve();
            return;
          }
          const clearedDF = removeEmptyCols(df);
          const replacer = this.valueReplacers.find((replacer) => replacer.check(tableName, colName, id));
          if (replacer) nameHost.textContent = replacer.replacer(clearedDF, id);

          ui.tooltip.bind(nameHost, () => ui.wait(async () => await this.renderDataFrame(clearedDF, schemaName, tableName)));
          nameHost.addEventListener('click', (e) => {
            e.stopImmediatePropagation();
            this.schemaInfoPromise().then((schemaAndConnection) => {
              if (schemaAndConnection && schemaAndConnection.connection && schemaAndConnection.schema)
                grok.shell.setCurrentObject(new DBValueObject(schemaAndConnection.connection.nqName, this.schemaName, tableName, colName, id, undefined, this.catalogName), false, true);
            });
          });
          nameHost.style.color = 'var(--blue-1)';
          nameHost.style.cursor = 'pointer';
          // after the successful load, try to add meaningful names to the meaningless ids
          const replaceCol = this.defaultHeaderReplacerColumns.find((col) => clearedDF.col(col) && !clearedDF.col(col)!.isNone(0));
          if (replaceCol) {
            // add id (replacer)
            nameHost.textContent = `${id} (${clearedDF.col(replaceCol)!.getString(0)})`;
          }
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

  async getTable(schemaName: string, tableName: string, match: string, matchValue: string | number, joinOptions: QueryJoinOptions[] = []) {
    const schemaAndConnection = await this.schemaInfoPromise();
    const res = await queryDB(schemaAndConnection?.connection ?? null, tableName, match, matchValue, schemaName, joinOptions, false, this.catalogName);
    return res;
  }

  async renderMultiRowTable(
    schemaName: string,
    tableName: string,
    match: string,
    matchValue: string | number,
    joinOptions: QueryJoinOptions[] = []
  ) {
    const df = await this.getTable(schemaName, tableName, match, matchValue, joinOptions);
    if (df.rowCount < 1)
      return {root: ui.divText('ID not found'), rowCount: 0};
    if (df.rowCount === 1)
      return {root: this.renderDataFrame(df, schemaName, tableName), rowCount: 1};
    const acc = ui.accordion(`Multiple rows for ${tableName}`);
    const dfMaxCount = Math.min(df.rowCount, 10);
    const replaceColName = ciGet(this.headerReplacers, tableName);
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
        return ui.wait(async () => await this.renderDataFrame(rowDf, schemaName, tableName));
      });
    }
    return {root: acc.root, rowCount: df.rowCount};
  }

  async renderTable(
    schemaName: string,
    tableName: string,
    match: string,
    matchValue: string | number,
    joinOptions: QueryJoinOptions[] = []
  ) {
    const df = await this.getTable(schemaName, tableName, match, matchValue, joinOptions);
    if (df.rowCount < 1) return ui.divText('ID not found');
    return this.renderDataFrame(df, schemaName, tableName);
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
  async renderDataFrame(df: DG.DataFrame, schemaName: string, tableName: string, options?: {keepEmptyValues: boolean, skipCustomSelected?: boolean}) {
    const schemaAndConnection = await this.schemaInfoPromise();
    if (!schemaAndConnection)
      return ui.divText('Schema information is not available');
    if (df.rowCount === 0)
      return ui.divText('No data');

    const clearedDF = options?.keepEmptyValues ? df: removeEmptyCols(df);
    await clearedDF.meta.detectSemanticTypes();
    let entries = clearedDF.columns.names()
      .map((colName) => [colName, clearedDF.col(colName)!.isNone(0) ? '' : clearedDF.col(colName)!.getString(0)]) as [string, string | number][];

    // reverse compatibility: first check schema qualified table name
    const customSelectedColumnEntries = ciGet(this.customSelectedColumns, `${schemaName}.${tableName}`) ?? ciGet(this.customSelectedColumns, tableName);
    if (!options?.skipCustomSelected && customSelectedColumnEntries) {
      const ciHas = (set: Set<string>, val: string) => {
        const lower = val.toLowerCase();
        for (const v of set) if (v.toLowerCase() === lower) return true;
        return false;
      };
      entries = entries.filter(([key, _]) => ciHas(customSelectedColumnEntries, key));
      // reorder entries according to customSelectedColumns if present
      const cols = Array.from(customSelectedColumnEntries);
      const colsLower = cols.map((c) => c.toLowerCase());
      entries.sort((a, b) => {
        return colsLower.indexOf(a[0]?.toLowerCase()) - colsLower.indexOf(b[0]?.toLowerCase());
      });
    }
    // similarly here, the column name in the df might not match the db column name if the df is a result of a query joining multiple tables
    const dbColumnNames: {[colName: string]: string} = {};
    const dbColTableNames: {[colName: string]: string} = {};
    const dbColSchemaNames: {[colName: string]: string} = {};
    entries.forEach(([colName, _]) => {
      const col = clearedDF.col(colName);
      if (col?.getTag(DG.Tags.DbColumn))
        dbColumnNames[colName] = col.getTag(DG.Tags.DbColumn)!;
      if (col?.getTag(DG.Tags.DbTable))
        dbColTableNames[colName] = col.getTag(DG.Tags.DbTable)!;
      if (col?.getTag(DG.Tags.DbSchema))
        dbColSchemaNames[colName] = col.getTag(DG.Tags.DbSchema)!;
    });

    const isAggregatedFromMultipleTables = Object.keys(dbColTableNames).some((colName) => !ciEquals(dbColTableNames[colName], tableName));

    const mainTable = ui.table(entries, (entry) => {
      const colName = entry[0];
      const dbColName = dbColumnNames[colName] ?? colName;
      const value = entry[1];
      const dbTableName = dbColTableNames[colName] ?? tableName; // here, it is mapped from the df column, not its db name
      const dbSchemaName = dbColSchemaNames[colName] ?? schemaName;
      const cutomRenderer = this.customRenderers.find((renderer) => renderer.check(dbTableName, dbColName, value));
      if (cutomRenderer) return [textRenderer(colName), cutomRenderer.renderer(value, schemaAndConnection.connection)];
      // if its autodetected as a molecule or helm or something else in the future
      const semTypeRenderer = this.semtypeRenderers.find((renderer) => {
        const col = clearedDF.col(colName);
        return renderer.check(col ?? undefined);
      });
      if (semTypeRenderer) return [textRenderer(colName), semTypeRenderer.renderer(value)];

      const refInfo = schemaAndConnection.schema.references?.[dbSchemaName]?.[dbTableName]?.[dbColName];
      if (refInfo)
        return [textRenderer(colName), this.refIdRenderer(schemaAndConnection.connection, value, refInfo.refSchema, refInfo.refTable, refInfo.refColumn)];
      // if the given column is referenced by other tables, then also use refIdRenderer
      const refedByInfo = schemaAndConnection.schema.referencedBy?.[dbSchemaName]?.[dbTableName]?.[dbColName];
      if (refedByInfo && refedByInfo.length > 0)
        return [textRenderer(colName), this.refIdRenderer(schemaAndConnection.connection, value, dbSchemaName, dbTableName, dbColName)];

      const uniqueColForTable = ciGet(this.uniqueColNames, dbTableName);
      const isUnqueCol = uniqueColForTable !== undefined && ciEquals(uniqueColForTable, dbColName);
      if (isUnqueCol) {
        if (!isAggregatedFromMultipleTables) // if its not aggregated, the tooltip in of the id will be the exact same as card, so we use simpler ownerIdRenderer
          return [textRenderer(colName), ownIdRenderer(value, dbTableName, dbColName, schemaAndConnection.connection.nqName, dbSchemaName, this.catalogName)];
        return [textRenderer(colName), this.refIdRenderer(schemaAndConnection.connection, value, dbSchemaName, dbTableName, dbColName)]; // otherwise use refIdRenderer
      }


      return [textRenderer(colName), textRenderer(value)];
    });
    return ui.divV([mainTable]);
  }

  async renderAssociations(acc: DG.Accordion, schemaPromise: () => Promise<SchemaAndConnection | null>, curSchema: string, curTable: string, curDf: DG.DataFrame) {
    const schemaAndConnection = await schemaPromise();
    if (!schemaAndConnection || (curDf?.rowCount ?? 0) === 0)
      return;
    this.names.useSchema(schemaAndConnection.schema);
    const groups = linkGroups(schemaAndConnection.schema, curSchema, curTable, curDf, this.entryPointOptions.joinOptions)
      .map((g) => ({...g, links: g.links.filter((l) => !this.names.entity(l.schema, l.table)?.hidden)}))
      .filter((g) => g.links.length > 0);
    if (groups.length === 0)
      return;
    const card: CardInfo = {schema: curSchema, table: curTable, df: curDf, connection: schemaAndConnection.connection};
    const linksPane = acc.addPane('Links', () => this.renderLinks(groups, card));
    const linksHeader = paneHeader(linksPane);
    if (linksHeader)
      ui.tooltip.bind(linksHeader, 'Rows of other tables that point at this entry, or at what it points at');
  }

  /** The Links section: a box for each key the rows point at, headed by what they share with the card, with one pane
   * for each table they are in */
  private renderLinks(groups: LinkGroup[], card: CardInfo): HTMLElement {
    const root = ui.divV([], 'db-explorer-links');
    Object.assign(root.style, {gap: '8px', width: '100%', minWidth: '0'});
    for (const g of groups) {
      const panes = this.linkPanes(g, card);
      // a wide card opened in a box scrolls inside it, and the box keeps the panel's width
      Object.assign(panes.style, {padding: '2px 4px 4px 14px', overflowX: 'auto'});
      const box = ui.divV([this.linkHeading(g, groups, card), panes], 'db-explorer-links-group');
      Object.assign(box.style, {border: '1px solid var(--grey-2)', borderRadius: '4px', maxWidth: '100%', minWidth: '0',
        boxSizing: 'border-box'});
      root.appendChild(box);
    }
    return root;
  }

  /** `This drug`, or `Same target · GROKTGT-000010`: what the rows of the group share with the card */
  private linkHeading(g: LinkGroup, groups: LinkGroup[], card: CardInfo): HTMLElement {
    const name = this.names.name(g.target.schema, g.target.table);
    let text: string;
    if (g.own)
      text = `This ${lcFirst(name.singular)}`;
    else if (g.label)
      text = `Same ${lcFirst(g.label)}`;
    else {
      const sameTable = (t: LinkRef) => ciEquals(t.table, g.target.table) && ciEquals(t.schema, g.target.schema);
      const twice = groups.filter((x) => !x.own && sameTable(x.target)).length > 1;
      const self = ciEquals(card.table, g.target.table) && ciEquals(card.schema, g.target.schema);
      const role = self || twice ? roleWords(g.field, this.names.words(g.target.schema, g.target.table, g.target.column)) : [];
      text = `Same ${role.length > 0 ? `${role.join(' ')} ` : ''}${lcFirst(name.singular)}`;
    }
    const label = ui.span([text]);
    ui.tooltip.bind(label, () => this.physicalTooltip(g.target.schema, g.target.table, g.target.column, g.value));
    const parts: HTMLElement[] = [label];
    if (!g.own) {
      const value = ui.span([String(g.value)]);
      Object.assign(value.style, {color: 'var(--blue-1)', cursor: 'pointer'});
      ui.tooltip.bind(value, () => ui.wait(async () => await this.renderTable(g.target.schema, g.target.table, g.target.column, g.value, this.entryPointOptions.joinOptions)));
      value.addEventListener('click', (e) => {
        e.stopImmediatePropagation();
        grok.shell.setCurrentObject(new DBValueObject(card.connection.nqName, g.target.schema, g.target.table, g.target.column, g.value, undefined, this.catalogName), false, true);
      });
      parts.push(ui.span([' · ']), value);
    }
    const heading = ui.divH(parts, 'db-explorer-links-heading');
    Object.assign(heading.style, {padding: '6px 8px', fontWeight: '500', color: 'var(--grey-6)', gap: '4px',
      alignItems: 'baseline', flexWrap: 'wrap', background: 'var(--grey-1)', borderBottom: '1px solid var(--grey-2)',
      borderRadius: '4px 4px 0 0'});
    heading.dataset.dbTarget = `${g.target.schema}.${g.target.table}.${g.target.column}`;
    return heading;
  }

  /** The table panes of a group: details first (rows that share the key), then records that have identifiers, then
   * results; views that repeat other tables are folded under More */
  private linkPanes(g: LinkGroup, card: CardInfo): HTMLElement {
    const acc = ui.accordion(`db-explorer-links:${g.target.schema}.${g.target.table}.${g.target.column}`);
    const items = g.links.map((l) => ({link: l, rank: this.linkRank(l), name: this.linkName(l, g)}));
    items.filter((it) => it.rank < LINK_RANK_MORE)
      .sort((a, b) => a.rank - b.rank)
      .forEach((it) => this.addLinkPane(acc, it.link, it.name, g, card));
    const more = items.filter((it) => it.rank === LINK_RANK_MORE);
    if (more.length > 0) {
      const morePane = acc.addPane('More', () => {
        const moreAcc = ui.accordion(`db-explorer-links-more:${g.target.schema}.${g.target.table}.${g.target.column}`);
        more.forEach((it) => this.addLinkPane(moreAcc, it.link, it.name, g, card));
        return moreAcc.root;
      });
      const header = paneHeader(morePane);
      if (header)
        ui.tooltip.bind(header, `Tables that repeat what the others hold: ${more.map((it) => it.name).join(', ')}`);
    }
    return acc.root;
  }

  /** 0: details that share the key, 1: records with a key of their own, 2: other rows (results), 3: views folded under More */
  private linkRank(l: LinkRef): number {
    const unique = ciGet(this.uniqueColNames, l.table);
    const known = {unique: !!unique, card: !!this.cardColumns(l.schema, l.table), header: !!ciGet(this.headerReplacers, l.table)};
    if (this.names.isMore(l.schema, l.table, known, Object.keys(this.customSelectedColumns).length))
      return LINK_RANK_MORE;
    const ep = this.names.entryPoint(l.table);
    if (!ep && (ciEquals(l.column, 'id') || (!!unique && ciEquals(l.column, unique))))
      return 0;
    return unique || ep ? 1 : 2;
  }

  /** `Assay results`, `Hierarchy · as parent`, `Child compounds (CPD)` */
  linkName(l: LinkRef, g: LinkGroup): string {
    if (l.reverseLabel)
      return l.reverseLabel;
    const n = this.names.name(l.schema, l.table);
    const targetWords = this.names.words(g.target.schema, g.target.table, g.target.column);
    let text = n.source === 'table' ? withoutParent(n.plural, targetWords) : n.plural;
    const sameTable = (x: LinkRef) => ciEquals(x.table, l.table) && ciEquals(x.schema, l.schema);
    const self = ciEquals(l.table, g.target.table) && ciEquals(l.schema, g.target.schema);
    const twice = g.links.filter(sameTable).length > 1;
    let suffix = '';
    if (self || twice) {
      const role = roleWords(l.column, targetWords);
      if (self && role.length === 1 && role[0] === 'parent')
        text = `Child ${lcFirst(text)}`;
      else if (self && role.length === 1 && role[0] === 'child')
        text = `Parent ${lcFirst(text)}`;
      else if (role.length > 0)
        suffix = ` · as ${role.join(' ')}`;
    }
    return `${text}${n.code ? ` (${n.code})` : ''}${suffix}`;
  }

  private physicalTooltip(schema: string, table: string, column: string, value: string | number): HTMLElement {
    const lines = [ui.divText(`${schema}.${table}.${column} = '${value}'`, {style: {fontFamily: 'var(--font-monospace, monospace)'}})];
    const description = this.names.description(schema, table);
    if (description)
      lines.push(ui.divText(description));
    return ui.divV(lines);
  }

  /** A pane of the rows of one table that point at the group's key. Its rows are read when it is opened, and their
   * count is shown on its header then. */
  private addLinkPane(acc: DG.Accordion, link: LinkRef, name: string, g: LinkGroup, card: CardInfo) {
    const count = ui.span([], 'db-explorer-links-count');
    Object.assign(count.style, {color: 'var(--grey-4)', marginLeft: '6px', fontWeight: 'normal', fontVariantNumeric: 'tabular-nums'});
    const addAll = () => this.addAllAssociated(link, g.value, `${name} of ${g.value}`);
    const pane = acc.addPane(name, () => ui.wait(async () => {
      const rows = await this.getTable(link.schema, link.table, link.column, g.value, this.entryPointOptions.joinOptions);
      count.textContent = rows.rowCount.toLocaleString('en-US');
      return this.renderLinkRows(rows, link, card, addAll);
    }));
    const header = paneHeader(pane);
    if (!header)
      return;
    header.dataset.dbLink = `${link.schema}.${link.table}.${link.column}`;
    header.appendChild(count);
    ui.tooltip.bind(header, () => this.physicalTooltip(link.schema, link.table, link.column, g.value));
    const title = `Add all ${name} to the workspace`;
    const icon = ui.icons.add(() => {}, title);
    // its own listener, so that the click does not also open the pane
    icon.addEventListener('click', (e) => {
      e.stopImmediatePropagation();
      e.preventDefault();
      addAll();
    });
    header.appendChild(icon);
    const menu = DG.Menu.popup().item(title, addAll);
    pane.root.addEventListener('contextmenu', (e) => {
      e.preventDefault();
      setTimeout(() => menu.show());
    });
  }

  private async addAllAssociated(link: LinkRef, value: string | number, viewName: string) {
    const pi = DG.TaskBarProgressIndicator.create('Opening all associated entries');
    try {
      const assocDf = await this.getTable(link.schema, link.table, link.column, value, this.entryPointOptions.joinOptions);
      if (assocDf.rowCount === 0) {
        grok.shell.info(`No associated ${link.table} entries found`);
        return;
      }
      assocDf.name = viewName;
      grok.shell.addTableView(assocDf);
    } catch (e) {
      console.error(e);
      grok.shell.error(`Failed to open associated ${link.table} entries`);
    } finally {
      pi.close();
    }
  }

  private cardColumns(schema: string, table: string): Set<string> | undefined {
    return ciGet(this.customSelectedColumns, `${schema}.${table}`) ?? ciGet(this.customSelectedColumns, table);
  }

  /** The table's key among the rows: its unique column, or the column of its identifier */
  private keyColumn(table: string, rows: DG.DataFrame): string | null {
    const name = ciGet(this.uniqueColNames, table) ?? this.names.entryPoint(table)?.column;
    return (name ? rows.col(name)?.name : null) ?? null;
  }

  /** The column that names a row in a list: the configured header column, else a column named like a name, else one
   * of the default header columns; [skip] are the key and the column every row shares */
  private titleColumn(table: string, rows: DG.DataFrame, skip: string[]): string | null {
    const usable = (name?: string | null): string | null => {
      const col = name ? rows.col(name) : null;
      if (!col || skip.some((s) => ciEquals(s, col.name)))
        return null;
      for (let i = 0; i < rows.rowCount; i++) {
        if (!col.isNone(i))
          return col.name;
      }
      return null;
    };
    return usable(ciGet(this.headerReplacers, table)) ??
      rows.columns.names().map((n) => NAME_COLUMN.test(n) ? usable(n) : null).find((n) => !!n) ??
      this.defaultHeaderReplacerColumns.map(usable).find((n) => !!n) ?? null;
  }

  /** The rows of a link, one line each: its key (a link to its card) and its name; the rest of the row opens under it */
  private async renderLinkRows(rows: DG.DataFrame, link: LinkRef, card: CardInfo, addAll: () => void): Promise<HTMLElement> {
    if (rows.rowCount === 0)
      return ui.divText('No records', {style: {color: 'var(--grey-4)'}});
    const shown = Math.min(rows.rowCount, MAX_MULTIROW_VALUES);
    const key = this.keyColumn(link.table, rows);
    // the column the rows point with holds the card's value in every row: it names none of them
    const skip = [link.column, ...(key ? [key] : [])];
    const titleCol = this.titleColumn(link.table, rows, skip) ??
      (key ? null : rows.columns.names().find((n) => !skip.some((s) => ciEquals(s, n)) && !rows.col(n)!.isNone(0)) ?? null);
    const keys = key ? Array.from({length: rows.rowCount}, (_, i) => rows.col(key)!.getString(i)) : [];
    const unique = key !== null && new Set(keys).size === keys.length;
    const cardKey = unique && ciEquals(link.table, card.table) && ciEquals(link.schema, card.schema) ? card.df.col(key!)?.getString(0) : undefined;

    const acc = ui.accordion();
    for (let i = 0; i < shown; i++) {
      const pane = acc.addPane('', () => {
        const bs = DG.BitSet.create(rows.rowCount);
        bs.set(i, true);
        return ui.wait(async () => await this.renderDataFrame(rows.clone(bs), link.schema, link.table));
      });
      const header = paneHeader(pane);
      if (!header)
        continue;
      const keyValue = key ? keys[i] : '';
      const title = titleCol && !rows.col(titleCol)!.isNone(i) ? rows.col(titleCol)!.getString(i) : '';
      header.dataset.dbRow = keyValue || title || `Row ${i + 1}`;
      const parts: HTMLElement[] = [];
      if (keyValue) {
        const id = ownIdRenderer(keyValue, link.table, key!, card.connection.nqName, link.schema, this.catalogName);
        // the key stays whole; a long name gives way
        Object.assign(id.style, {flexShrink: '0', maxWidth: 'none'});
        parts.push(id);
      }
      if (title && title !== keyValue) {
        const t = textRenderer(title, false);
        Object.assign(t.style, {flexShrink: '1', minWidth: '0'});
        if (keyValue)
          t.style.color = 'var(--grey-4)';
        parts.push(t);
      }
      if (parts.length === 0)
        parts.push(ui.span([`Row ${i + 1}`]));
      if (cardKey !== undefined && keyValue === cardKey)
        parts.push(thisTag());
      const line = ui.divH(parts);
      Object.assign(line.style, {gap: '8px', alignItems: 'baseline', display: 'inline-flex', minWidth: '0'});
      header.appendChild(line);
    }
    const root = ui.divV([acc.root], 'db-explorer-link-rows');
    if (rows.rowCount > shown) {
      const openAll = ui.span(['Open all as a table']);
      Object.assign(openAll.style, {color: 'var(--blue-1)', cursor: 'pointer'});
      openAll.addEventListener('click', addAll);
      const footer = ui.divH([ui.span([`Showing ${shown} of ${rows.rowCount.toLocaleString('en-US')}`]), ui.span(['·']), openAll]);
      Object.assign(footer.style, {color: 'var(--grey-4)', fontSize: '12px', paddingTop: '4px', gap: '4px'});
      root.appendChild(footer);
    }
    return root;
  }
}

type CardInfo = {schema: string, table: string, df: DG.DataFrame, connection: DG.DataConnection};

const LINK_RANK_MORE = 3;
/** A column that names a row: `name`, `alias`, `pref_name`, `compound_name`, `ALIAS` */
const NAME_COLUMN = /^(name|alias|title|label|pref_name|display_name)$|_(name|alias)$/i;

function paneHeader(pane: DG.AccordionPane): HTMLElement | null {
  return pane.root.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
}

function thisTag(): HTMLElement {
  const tag = ui.span(['this'], 'db-explorer-this');
  Object.assign(tag.style, {fontSize: '10px', border: '1px solid var(--grey-2)', borderRadius: '3px', padding: '0 4px', color: 'var(--grey-4)'});
  ui.tooltip.bind(tag, 'The entry this card is about');
  return tag;
}


export class ExampleExplorerRenderer extends DBExplorerRenderer {
  private static schemaCache: Map<string, Promise<SchemaAndConnection | null>> = new Map();

  constructor(
    connection: DG.DataConnection,
    shemaName: string,
    entryPointOptions: EntryPointOptions,
    catalogName?: string
  ) {
    const cacheKey = `${connection.id}||${catalogName ?? ''}||${shemaName}`;
    let schemaPromise = ExampleExplorerRenderer.schemaCache.get(cacheKey);
    if (!schemaPromise) {
      schemaPromise = getSchemaReferences(connection, shemaName, catalogName).then((schema) => ({schema, connection}));
      ExampleExplorerRenderer.schemaCache.set(cacheKey, schemaPromise);
    }
    super(() => schemaPromise!, shemaName, entryPointOptions, catalogName);
  }

  override async renderDataFrame(df: DG.DataFrame, schemaName: string, tableName: string, _opts?: {keepEmptyValues: boolean, columnTablesMap?: {[colName: string]: string}, skipCustomSelected?: boolean}): Promise<HTMLDivElement> {
    return super.renderDataFrame(df, schemaName, tableName, {keepEmptyValues: true});
  }

  override async renderTable(schemaName: string, tableName: string, match: string, matchValue: string | number, joinOptions?: QueryJoinOptions[]): Promise<HTMLDivElement> {
    const schemaAndConnection = await this.schemaInfoPromise();
    if (!schemaAndConnection)
      return ui.divText('Schema information is not available');
    const df = await queryDB(schemaAndConnection.connection, tableName, '', '', this.schemaName, joinOptions, true, this.catalogName);
    if (df.rowCount < 1) return ui.divText('ID not found');
    return this.renderDataFrame(df, schemaName, tableName);
  }
}

export function renderExampleCard(connection: DG.DataConnection, schemaName: string, tableName: string,
  config: Partial<DBExplorerConfig>): HTMLElement {
  const ex = new ExampleExplorerRenderer(
    connection,
    schemaName,
    {joinOptions: config.joinOptions ?? [], valueConverter: (a) => a},
    config.catalogName
  );
  ex.addCustomSelectedColumns(config.customSelectedColumns ?? {});
  (config.customRenderers ?? []).forEach((cr) =>
    ex.addCustomRenderer((t, c, _v) => t === cr.table && c === cr.column, (v) => getDefaultRendererByName(cr.renderer)(v))
  );
  (config.uniqueColumns) && ex.addUniqueColNames(config.uniqueColumns);
  (config.customSelectedColumns) && ex.addCustomSelectedColumns(config.customSelectedColumns);
  const c = ui.card(ui.wait(async () => {
    return await ex.renderTable(schemaName, tableName, '', '', config.joinOptions ?? []);
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
