import JSZip from 'jszip';
import {MwxWorksheet} from './mwx-types';
import {parseSheetColumns} from '../shared/sheet-parser';


/** Parses MWX (Minitab Worksheet) bytes into a typed worksheet model. */
export async function parseMwx(bytes: Uint8Array): Promise<MwxWorksheet> {
  const zip = await JSZip.loadAsync(bytes);

  const metaFile = zip.file('sheet_metadata_20.json') ?? zip.file('/sheet_metadata_20.json');
  let wsName = 'Worksheet';
  let version = 0;
  if (metaFile) {
    const meta = JSON.parse(await metaFile.async('text'));
    wsName = meta.Worksheet?.Name ?? wsName;
    version = meta.Version ?? 0;
  }

  const sheetFile = zip.file('sheets/0/sheet.json') ?? zip.file('/sheets/0/sheet.json');
  if (!sheetFile)
    throw new Error('MWX file does not contain sheets/0/sheet.json');

  const sheet = JSON.parse(await sheetFile.async('text'));
  const data = sheet.Data ?? sheet;
  const title = data.PrivateTitle ?? wsName;
  const {columns, rowCount} = parseSheetColumns(data);

  return {name: title, version, columns, rowCount};
}
