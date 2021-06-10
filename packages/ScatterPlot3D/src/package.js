/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
//import {test222} from "./test.js"
import {ScatterPlot3Dviewer} from './ScatterPlot3Dviewer.js';
 
export let _package = new DG.Package();

//name: test
export function test() {
  grok.shell.info(_package.webRoot);
}

console.error('bbbb11') 
 
//name: ScatterPlot3D
//description: Creates an awesome viewer
//tags: viewer
//output: viewer result
export function scatterPlot3Dviewer() {
  return new ScatterPlot3Dviewer();
}
