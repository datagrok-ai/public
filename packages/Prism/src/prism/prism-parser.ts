import JSZip from 'jszip';

import {
  PrismFile, PrismSheet, PrismAnalysis,
  DocumentJson, SheetJson, DatasetJson,
} from './prism-types';


/** Reads and parses a JSON file from the ZIP archive by exact path. */
async function readJsonAt<T>(zip: JSZip, path: string): Promise<T | null> {
  const file = zip.file(path);
  if (!file)
    return null;
  const text = await file.async('string');
  return JSON.parse(text) as T;
}

/** Finds a file by path pattern in the ZIP. */
function findFile(zip: JSZip, pattern: RegExp): string | null {
  for (const path of Object.keys(zip.files)) {
    if (pattern.test(path))
      return path;
  }
  return null;
}

/** Parses a CSV string into a 2D array, handling quoted fields. */
function parseCsv(text: string): string[][] {
  const rows: string[][] = [];
  const lines = text.split('\n');
  for (const line of lines) {
    const trimmed = line.replace(/\r$/, '');
    if (trimmed === '')
      continue;

    // Handle quoted fields
    const fields: string[] = [];
    let i = 0;
    while (i <= trimmed.length) {
      if (i === trimmed.length) {
        fields.push('');
        break;
      }
      if (trimmed[i] === '"') {
        let j = i + 1;
        let val = '';
        while (j < trimmed.length) {
          if (trimmed[j] === '"' && j + 1 < trimmed.length && trimmed[j + 1] === '"') {
            val += '"';
            j += 2;
          }
          else if (trimmed[j] === '"') {
            j++;
            break;
          }
          else {
            val += trimmed[j];
            j++;
          }
        }
        fields.push(val);
        i = j + 1; // skip comma
      }
      else {
        const commaIdx = trimmed.indexOf(',', i);
        if (commaIdx === -1) {
          fields.push(trimmed.substring(i));
          break;
        }
        fields.push(trimmed.substring(i, commaIdx));
        i = commaIdx + 1;
      }
    }
    rows.push(fields);
  }
  return rows;
}


/** Reads the title from a dataset JSON at data/sets/<datasetId>.json. */
async function readDatasetTitle(zip: JSZip, datasetId: string): Promise<string> {
  const json = await readJsonAt<DatasetJson>(zip, `data/sets/${datasetId}.json`);
  return json?.title ?? '';
}

/** Reads a sheet.json for a data sheet at data/sheets/<sheetId>/sheet.json. */
async function readDataSheetJson(zip: JSZip, sheetId: string): Promise<SheetJson | null> {
  return readJsonAt<SheetJson>(zip, `data/sheets/${sheetId}/sheet.json`);
}

/** Reads a sheet.json for an analysis at analyses/<sheetId>/sheet.json. */
async function readAnalysisSheetJson(zip: JSZip, sheetId: string): Promise<any | null> {
  return readJsonAt<any>(zip, `analyses/${sheetId}/sheet.json`);
}


/** Parses a data sheet: reads sheet.json, data.csv, and dataset metadata. */
async function parseDataSheet(zip: JSZip, sheetId: string): Promise<PrismSheet | null> {
  const sheetJson = await readDataSheetJson(zip, sheetId);
  if (!sheetJson?.table)
    return null;

  const tableUid = sheetJson.table.uid;

  // Read CSV from data/tables/<tableUid>/data.csv
  const csvPath = `data/tables/${tableUid}/data.csv`;
  const csvFile = zip.file(csvPath);
  if (!csvFile)
    return null;

  const csvText = await csvFile.async('string');
  const data = parseCsv(csvText);

  // Read Y column titles from data/sets/<datasetId>.json
  const yColumnTitles: string[] = [];
  if (sheetJson.table.dataSets) {
    for (const dsId of sheetJson.table.dataSets) {
      const title = await readDatasetTitle(zip, dsId);
      yColumnTitles.push(title);
    }
  }

  // Read X column title
  let xTitle: string | null = null;
  if (sheetJson.table.xDataSet)
    xTitle = await readDatasetTitle(zip, sheetJson.table.xDataSet) || null;

  return {
    id: sheetId,
    title: sheetJson.title ?? '',
    dataFormat: sheetJson.table.dataFormat ?? 'y_single',
    replicatesCount: sheetJson.table.replicatesCount ?? 1,
    xTitle,
    rowTitlesPresent: !!sheetJson.table.rowTitlesDataSet,
    yColumnTitles,
    data,
  };
}


/** Parses an analysis sheet: follows preferredResultSheet chain to get result data. */
async function parseAnalysisSheet(
  zip: JSZip, sheetId: string, sheetName: string
): Promise<PrismAnalysis | null> {
  const sheetJson = await readAnalysisSheetJson(zip, sheetId);
  if (!sheetJson?.preferredResultSheet)
    return null;

  // Read the result sheet reference from analyses/<id>/result_sheets/<resultId>.json
  const resultId = sheetJson.preferredResultSheet;
  const resultRefPath = `analyses/${sheetId}/result_sheets/${resultId}.json`;
  const resultRef = await readJsonAt<any>(zip, resultRefPath);
  if (!resultRef?.dataSheet)
    return null;

  // Read the actual data sheet for the analysis results
  const resultSheetJson = await readDataSheetJson(zip, resultRef.dataSheet);
  if (!resultSheetJson?.table)
    return null;

  const csvPath = `data/tables/${resultSheetJson.table.uid}/data.csv`;
  const csvFile = zip.file(csvPath);
  if (!csvFile)
    return null;

  const csvText = await csvFile.async('string');
  const resultData = parseCsv(csvText);

  // Read column titles from datasets
  const columnTitles: string[] = [];
  if (resultSheetJson.table.dataSets) {
    for (const dsId of resultSheetJson.table.dataSets) {
      const title = await readDatasetTitle(zip, dsId);
      columnTitles.push(title);
    }
  }

  return {
    id: sheetId,
    title: sheetName,
    resultData,
    columnTitles,
  };
}


/** Extracts sheet name from sheetAttributesMap entry (value is {title: "..."} or string). */
function extractSheetName(value: any): string {
  if (typeof value === 'string')
    return value;
  if (value && typeof value === 'object' && value.title)
    return value.title;
  return '';
}


/** Parses a .prism file (ZIP archive) and returns a structured PrismFile object. */
export async function parsePrismFile(bytes: Uint8Array): Promise<PrismFile> {
  const zip = await JSZip.loadAsync(bytes);

  // Read top-level document.json
  const doc = await readJsonAt<DocumentJson>(zip, 'document.json');
  if (!doc)
    return {sheets: [], analyses: [], sheetNames: new Map()};

  // Build sheet name map (values can be {title: "..."} objects or plain strings)
  const sheetNames = new Map<string, string>();
  if (doc.sheetAttributesMap) {
    for (const [id, value] of Object.entries(doc.sheetAttributesMap))
      sheetNames.set(id, extractSheetName(value));
  }

  // Parse data sheets
  const sheets: PrismSheet[] = [];
  const dataSheetIds = doc.sheets?.data ?? [];
  for (const sheetId of dataSheetIds) {
    const sheet = await parseDataSheet(zip, sheetId);
    if (sheet) {
      if (!sheet.title && sheetNames.has(sheetId))
        sheet.title = sheetNames.get(sheetId)!;
      sheets.push(sheet);
    }
  }

  // Parse analysis sheets
  const analyses: PrismAnalysis[] = [];
  const analysisSheetIds = doc.sheets?.analyses ?? [];
  for (const sheetId of analysisSheetIds) {
    const name = sheetNames.get(sheetId) ?? sheetId;
    const analysis = await parseAnalysisSheet(zip, sheetId, name);
    if (analysis)
      analyses.push(analysis);
  }

  return {sheets, analyses, sheetNames};
}
