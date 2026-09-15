import {TableauPackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//description: Opens Tableau workbook file
//input: list bytes 
//output: list<dataframe> result
//meta.role: FileHandler
//meta.ext: twb
export function importTwb(bytes: Uint8Array) : any {
  return TableauPackageFunctions.importTwb(bytes);
}

//input: file file 
//output: view result
//meta.role: FileViewer
//meta.fileViewer: twb
export async function previewTwb(file: DG.FileInfo) : Promise<any> {
  return await TableauPackageFunctions.previewTwb(file);
}
