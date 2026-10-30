import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//meta.role: init
export async function init() : Promise<void> {
  await PackageFunctions.init();
}

//meta.autostartImmediate: true
//meta.role: autostart
export function emsAutostart() : void {
  PackageFunctions.emsAutostart();
}

//name: MolTrack
//description: Registers and searches compounds, batches and assay data
//input: string path { meta.url: true; optional: true; description: View to open: Compounds, Batches, Assays, AssayRuns, AssayResults or Additions }
//output: view result
//meta.browsePath: Chem
//meta.icon: images/moltrack.png
//meta.role: app
export async function molTrackApp(path: string) : Promise<any> {
  return await PackageFunctions.molTrackApp(path);
}

//input: dynamic appNode 
//meta.role: appTreeBrowser
//meta.app: MolTrack
export async function molTrackAppTreeBrowser(appNode: any) : Promise<void> {
  await PackageFunctions.molTrackAppTreeBrowser(appNode);
}

//description: Loads the mol-track demo set into the moltrack domain schema where its tables are empty
//output: string result
export async function loadEmsDemoData() : Promise<string> {
  return await PackageFunctions.loadEmsDemoData();
}
