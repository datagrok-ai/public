import {TwbFile, TwbDatasource, TwbColumn, TwbWorksheet} from './tableau-types';


function cleanName(name: string): string {
  return name.replace(/^\[/, '').replace(/\]$/, '');
}


function parseDatasources(root: Element): TwbDatasource[] {
  const results: TwbDatasource[] = [];
  const dsElements = root.querySelectorAll('datasources > datasource');

  for (const dsEl of Array.from(dsElements)) {
    const name = dsEl.getAttribute('name') || '';
    const caption = dsEl.getAttribute('caption') || '';

    // Extract source file info
    let sourceFile = '';
    let sourceDirectory = '';
    const connEl = dsEl.querySelector('connection > named-connections > named-connection > connection');
    if (connEl) {
      sourceFile = connEl.getAttribute('filename') || '';
      sourceDirectory = connEl.getAttribute('directory') || '';
    }

    // Build column map from <relation><columns><column>
    const columnMap = new Map<string, TwbColumn>();
    const relCols = dsEl.querySelectorAll('relation > columns > column');
    for (const col of Array.from(relCols)) {
      const colName = col.getAttribute('name') || '';
      columnMap.set(colName, {
        name: colName,
        caption: colName,
        datatype: col.getAttribute('datatype') || '',
        ordinal: parseInt(col.getAttribute('ordinal') || '0', 10),
        role: '',
        type: '',
        aggregation: '',
        containsNull: false,
      });
    }

    // Enrich from metadata-records
    const metaRecords = dsEl.querySelectorAll('metadata-records > metadata-record');
    for (const rec of Array.from(metaRecords)) {
      if (rec.getAttribute('class') !== 'column')
        continue;
      const remoteName = rec.querySelector('remote-name')?.textContent || '';
      const col = columnMap.get(remoteName);
      if (!col)
        continue;
      col.aggregation = rec.querySelector('aggregation')?.textContent || '';
      col.containsNull = rec.querySelector('contains-null')?.textContent === 'true';
    }

    // Enrich from top-level <column> children of datasource (caption, role, type)
    const summaryCols = Array.from(dsEl.children).filter(
      (el) => el.tagName === 'column' && el.getAttribute('name')
    );
    for (const el of summaryCols) {
      const rawName = cleanName(el.getAttribute('name') || '');
      if (rawName.includes('__tableau_internal'))
        continue;
      const col = columnMap.get(rawName);
      if (!col)
        continue;
      const cap = el.getAttribute('caption');
      if (cap)
        col.caption = cap;
      col.role = el.getAttribute('role') || col.role;
      col.type = el.getAttribute('type') || col.type;
    }

    // Filter out internal columns and sort by ordinal
    const columns = Array.from(columnMap.values())
      .filter((c) => !c.name.includes('__tableau_internal'))
      .sort((a, b) => a.ordinal - b.ordinal);

    results.push({name, caption, sourceFile, sourceDirectory, columns});
  }

  return results;
}


function parseWorksheets(root: Element): TwbWorksheet[] {
  const results: TwbWorksheet[] = [];
  const wsElements = root.querySelectorAll('worksheets > worksheet');

  for (const wsEl of Array.from(wsElements)) {
    const name = wsEl.getAttribute('name') || '';

    // Extract datasource name from datasource-dependencies
    let datasourceName = '';
    const depEl = wsEl.querySelector('datasource-dependencies');
    if (depEl)
      datasourceName = depEl.getAttribute('datasource') || '';

    // Collect used column names
    const usedColumns: string[] = [];
    const depCols = wsEl.querySelectorAll('datasource-dependencies > column');
    for (const col of Array.from(depCols)) {
      const colName = cleanName(col.getAttribute('name') || '');
      if (colName && !colName.includes('__tableau_internal'))
        usedColumns.push(colName);
    }

    // Read rows/cols expressions
    const rows = wsEl.querySelector('table > rows')?.textContent || '';
    const cols = wsEl.querySelector('table > cols')?.textContent || '';

    // Read mark class
    const markClass = wsEl.querySelector('pane > mark')?.getAttribute('class') || '';

    results.push({name, datasourceName, rows, cols, markClass, usedColumns});
  }

  return results;
}


export function parseTwbFile(text: string): TwbFile {
  const parser = new DOMParser();
  const doc = parser.parseFromString(text, 'text/xml');
  const root = doc.documentElement;

  return {
    version: root.getAttribute('version') || '',
    datasources: parseDatasources(root),
    worksheets: parseWorksheets(root),
  };
}
