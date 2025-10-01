import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: info
export function info() : void {
  PackageFunctions.info();
}

//name: Excalidraw
//tags: app
//output: view v
//meta.browsePath: Misc
//meta.icon: images/excal.png
export function excalidrawApp() {
  return PackageFunctions.excalidrawApp();
}

//description: Excalidraw viewer
//tags: fileViewer
//input: file fileContent 
//output: view result
//meta.ext: excalidraw
//meta.fileViewer: excalidraw
export async function excalfileViewer(fileContent: DG.FileInfo) : Promise<any> {
  return await PackageFunctions.excalfileViewer(fileContent);
}
