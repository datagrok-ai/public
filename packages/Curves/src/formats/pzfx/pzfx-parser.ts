import * as DG from 'datagrok-api/dg';

import {
  IFitChartData,
  IFitChartOptions,
  IFitSeries,
  IFitPoint,
  FIT_FUNCTION_SIGMOID,
} from '@datagrok-libraries/statistics/src/fit/fit-curve';


const PZFX_NS = 'http://graphpad.com/prism/Prism.htm';


export interface PzfxSubcolumn {
  values: (number | null)[];
  excluded: boolean[];
}

export interface PzfxYColumn {
  title: string;
  subcolumns: PzfxSubcolumn[];
}

export interface PzfxTable {
  id: string;
  title: string;
  tableType: string;
  xFormat: string;
  xTitle: string;
  xValues: (number | null)[];
  yColumns: PzfxYColumn[];
}


/** Returns child elements matching a local name, handling XML namespaces. */
function getChildren(parent: Element, localName: string): Element[] {
  const results: Element[] = [];
  for (let i = 0; i < parent.children.length; i++) {
    const child = parent.children[i];
    if (child.localName === localName)
      results.push(child);
  }
  return results;
}

/** Returns the first child element matching a local name. */
function getChild(parent: Element, localName: string): Element | null {
  for (let i = 0; i < parent.children.length; i++) {
    if (parent.children[i].localName === localName)
      return parent.children[i];
  }
  return null;
}

/** Returns the text content of the first child element matching a local name. */
function getChildText(parent: Element, localName: string): string {
  const el = getChild(parent, localName);
  return el?.textContent?.trim() ?? '';
}

/** Parses `<d>` elements from a `<Subcolumn>` element. */
function parseSubcolumn(subcolEl: Element): PzfxSubcolumn {
  const dElements = getChildren(subcolEl, 'd');
  const values: (number | null)[] = [];
  const excluded: boolean[] = [];

  for (const d of dElements) {
    const text = d.textContent?.trim() ?? '';
    if (text === '') {
      values.push(null);
      excluded.push(false);
    } else {
      values.push(parseFloat(text));
      excluded.push(d.getAttribute('Excluded') === '1');
    }
  }

  return {values, excluded};
}

/** Parses a `<YColumn>` element. */
function parseYColumn(yColEl: Element): PzfxYColumn {
  const title = getChildText(yColEl, 'Title');
  const subcolumns = getChildren(yColEl, 'Subcolumn').map(parseSubcolumn);
  return {title, subcolumns};
}

/** Parses a `<Table>` or `<HugeTable>` element. */
function parseTable(tableEl: Element): PzfxTable {
  const id = tableEl.getAttribute('ID') ?? '';
  const tableType = tableEl.getAttribute('TableType') ?? '';
  const xFormat = tableEl.getAttribute('XFormat') ?? '';
  const title = getChildText(tableEl, 'Title');

  let xTitle = '';
  let xValues: (number | null)[] = [];

  const xCol = getChild(tableEl, 'XColumn');
  if (xCol) {
    xTitle = getChildText(xCol, 'Title');
    const xSubcol = getChild(xCol, 'Subcolumn');
    if (xSubcol)
      xValues = parseSubcolumn(xSubcol).values;
  }

  const yColumns = getChildren(tableEl, 'YColumn').map(parseYColumn);

  return {id, title, tableType, xFormat, xTitle, xValues, yColumns};
}


/** Parses a PZFX XML string and returns an array of tables.
 * Truncates at the closing root tag to handle trailing binary data. */
export function parsePzfxXml(xmlText: string): PzfxTable[] {
  const closingTag = '</GraphPadPrismFile>';
  const closingIdx = xmlText.indexOf(closingTag);
  if (closingIdx !== -1)
    xmlText = xmlText.substring(0, closingIdx + closingTag.length);

  const parser = new DOMParser();
  const doc = parser.parseFromString(xmlText, 'text/xml');

  const root = doc.documentElement;
  if (root.nodeName === 'parsererror')
    return [];

  const tables: PzfxTable[] = [];
  for (let i = 0; i < root.children.length; i++) {
    const el = root.children[i];
    if (el.localName === 'Table' || el.localName === 'HugeTable')
      tables.push(parseTable(el));
  }
  return tables;
}


/** Converts a PZFX XY table to an IFitChartData object.
 * Returns null for non-XY tables or tables without numeric X data. */
export function pzfxTableToFitChartData(table: PzfxTable): IFitChartData | null {
  if (table.tableType !== 'XY' || table.xValues.length === 0)
    return null;

  const series: IFitSeries[] = [];

  for (const yCol of table.yColumns) {
    const points: IFitPoint[] = [];

    for (const subcol of yCol.subcolumns) {
      for (let i = 0; i < table.xValues.length; i++) {
        const x = table.xValues[i];
        if (x == null || i >= subcol.values.length)
          continue;
        const y = subcol.values[i];
        if (y == null)
          continue;
        points.push({
          x,
          y,
          outlier: subcol.excluded[i] || false,
        });
      }
    }

    if (points.length === 0)
      continue;

    series.push({
      name: yCol.title,
      fitFunction: FIT_FUNCTION_SIGMOID,
      showFitLine: true,
      showPoints: 'points',
      clickToToggle: true,
      droplines: ['IC50'],
      points,
    });
  }

  if (series.length === 0)
    return null;

  const chartOptions: IFitChartOptions = {
    logX: true,
    xAxisName: table.xTitle || 'Concentration',
    yAxisName: 'Response',
    title: table.title,
  };

  return {chartOptions, series};
}


/** Converts PZFX XY tables into a DataFrame with a fit column. Each row is one XY table. */
export function pzfxToFitDataFrame(tables: PzfxTable[]): DG.DataFrame {
  const xyTables = tables.filter((t) => t.tableType === 'XY');
  const names: string[] = [];
  const curves: string[] = [];

  for (const table of xyTables) {
    const chartData = pzfxTableToFitChartData(table);
    if (!chartData)
      continue;
    names.push(table.title || table.id);
    curves.push(JSON.stringify(chartData));
  }

  const df = DG.DataFrame.create(names.length);
  df.setTag(DG.Tags.Id, crypto.randomUUID());
  df.name = 'PZFX Curves';

  const nameCol = df.columns.addNewString('Table');
  for (let i = 0; i < names.length; i++)
    nameCol.set(i, names[i]);

  const curveCol = df.columns.addNewString('Fitted Curve');
  curveCol.semType = 'fit';
  curveCol.meta.cellRenderer = 'fit';
  for (let i = 0; i < curves.length; i++)
    curveCol.set(i, curves[i]);

  return df;
}


/** Converts a non-XY PZFX table into a plain DataFrame. */
export function pzfxTableToDataFrame(table: PzfxTable): DG.DataFrame {
  const maxRows = Math.max(...table.yColumns.map((yc) =>
    Math.max(...yc.subcolumns.map((sc) => sc.values.length), 0)), 0);

  const df = DG.DataFrame.create(maxRows);
  df.setTag(DG.Tags.Id, crypto.randomUUID());
  df.name = table.title || table.id;

  for (const yCol of table.yColumns) {
    if (yCol.subcolumns.length === 1) {
      const col = df.columns.addNewFloat(yCol.title || 'Values');
      const sc = yCol.subcolumns[0];
      for (let i = 0; i < sc.values.length; i++)
        col.set(i, sc.values[i] ?? DG.FLOAT_NULL);
    } else {
      for (let si = 0; si < yCol.subcolumns.length; si++) {
        const col = df.columns.addNewFloat(`${yCol.title}_${si + 1}`);
        const sc = yCol.subcolumns[si];
        for (let i = 0; i < sc.values.length; i++)
          col.set(i, sc.values[i] ?? DG.FLOAT_NULL);
      }
    }
  }

  return df;
}
