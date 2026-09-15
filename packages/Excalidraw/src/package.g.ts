import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: info
export function info() : void {
  PackageFunctions.info();
}

//name: Excalidraw
//output: view v
//meta.role: App
//meta.browsePath: Misc
//meta.icon: images/excal.png
export function excalidrawApp() {
  return PackageFunctions.excalidrawApp();
}

//description: Excalidraw viewer
//input: file fileContent 
//output: view result
//meta.ext: excalidraw
//meta.role: FileViewer
//meta.fileViewer: excalidraw
export async function excalfileViewer(fileContent: DG.FileInfo) : Promise<any> {
  return await PackageFunctions.excalfileViewer(fileContent);
}
