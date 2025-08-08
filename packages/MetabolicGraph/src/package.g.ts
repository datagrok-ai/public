import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: MetabolicGraph
//description: Metabolic graph application
//tags: app
//input: string path { meta.url: true; optional: true }
//input: string filter { optional: true }
//output: view result
//meta.icon: files/icons/metabolic.png
//meta.browsePath: Misc
export function metabolicGraphApp(path?: string, filter?: string) {
  return PackageFunctions.metabolicGraphApp(path, filter);
}

//name: EscherFileViewer
//tags: fileViewer
//input: file file 
//output: view result
//meta.fileViewer: json
//meta.fileViewerCheck: Metabolicgraph:escherFileViewerCheck
export async function escherFileViewer(file: DG.FileInfo) {
  return PackageFunctions.escherFileViewer(file);
}

//name: escherFileViewerCheck
//input: string content 
//output: bool result
export async function escherFileViewerCheck(content: string) {
  return PackageFunctions.escherFileViewerCheck(content);
}
