import JSZip from 'jszip';
import {MwxWorksheet} from '../mwx/mwx-types';
import {MpxProject} from './mpx-types';
import {parseSheetColumns} from '../shared/sheet-parser';


/** Parses MPX (Minitab Project) bytes into a typed project model. */
export async function parseMpx(bytes: Uint8Array): Promise<MpxProject> {
  const zip = await JSZip.loadAsync(bytes);

  // Read project metadata
  const metaFile = zip.file('project_metadata_20.json') ?? zip.file('/project_metadata_20.json');
  let meta: any = {};
  if (metaFile)
    meta = JSON.parse(await metaFile.async('text'));

  const creator = meta.Creator ?? '';
  const comments = meta.Comments ?? '';
  const history: string[] = meta.History ?? [];

  // Build worksheet name/ordering map from metadata
  const wsItems: any[] = meta.Worksheets?.Items ?? [];
  const wsOrdering: string[] = meta.Worksheets?.ExtraInformation?.worksheet_ordering ?? [];
  const wsNameById = new Map<string, string>();
  for (const item of wsItems)
    wsNameById.set(item.Id, item.Name ?? 'Worksheet');

  // Find all sheet.json files
  const sheetPaths = Object.keys(zip.files)
    .filter((f) => /^\/?sheets\/\d+\/sheet\.json$/.test(f))
    .sort();

  // Parse each worksheet
  const sheetsByIndex = new Map<number, MwxWorksheet>();
  const sheetIdByIndex = new Map<number, string>();

  for (const path of sheetPaths) {
    const match = path.match(/sheets\/(\d+)\/sheet\.json/);
    if (!match)
      continue;
    const index = parseInt(match[1]);
    const file = zip.file(path);
    if (!file)
      continue;

    const sheetJson = JSON.parse(await file.async('text'));
    const data = sheetJson.Data ?? sheetJson;
    const title = data.PrivateTitle ?? data.Name ?? `Worksheet ${index + 1}`;
    const wsId = data.WorksheetId_DEP ?? '';
    const {columns, rowCount} = parseSheetColumns(data);

    sheetsByIndex.set(index, {name: title, version: 0, columns, rowCount});
    if (wsId)
      sheetIdByIndex.set(index, wsId);
  }

  // Order worksheets: prefer metadata ordering, fall back to numeric index
  let worksheets: MwxWorksheet[];
  if (wsOrdering.length > 0) {
    const idToSheet = new Map<string, MwxWorksheet>();
    for (const [idx, ws] of sheetsByIndex) {
      const id = sheetIdByIndex.get(idx) ?? '';
      if (id)
        idToSheet.set(id, ws);
    }
    worksheets = [];
    for (const id of wsOrdering) {
      const ws = idToSheet.get(id);
      if (ws)
        worksheets.push(ws);
    }
    // Append any sheets not in the ordering
    for (const [, ws] of sheetsByIndex) {
      if (!worksheets.includes(ws))
        worksheets.push(ws);
    }
  }
  else {
    worksheets = [...sheetsByIndex.entries()]
      .sort((a, b) => a[0] - b[0])
      .map(([, ws]) => ws);
  }

  return {creator, comments, history, worksheets};
}
