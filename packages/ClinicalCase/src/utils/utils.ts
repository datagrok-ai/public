import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {scripts} from '../package-api';

export function updateDivInnerHTML(div: HTMLElement, content: any) {
  div.innerHTML = '';
  div.append(content);
}

export function createFilters(df: DG.DataFrame) {
  return DG.Viewer.fromType('Filters', df, {
    'showContextMenu': false,
  });
}

export function removeExtension(filename: string): string {
  const lastDotIndex = filename.lastIndexOf('.');
  return lastDotIndex === -1 ? filename : filename.substring(0, lastDotIndex);
};

export async function readClinicalFile(file: DG.FileInfo): Promise<DG.DataFrame> {
  let df: DG.DataFrame | null = null;
  try {
    if (file.extension === 'xpt')
      df = await scripts.readSas(file);
    else
      df = await DG.DataFrame.fromCsv(await file.readAsString());
  } catch (e: any) {
    grok.shell.error(`Error loading ${file.name}: ${e?.message ?? e}`);
  }
  return df;
}
