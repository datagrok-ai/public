/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { BoltzService } from './utils/boltz-service';
import { Boltz1AppView } from './demo/boltz-app';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: getBoltzConfigFolders
//output: list<string> configFiles
export async function getBoltzConfigFolders(): Promise<string[]> {
  return await BoltzService.getBoltzConfigFolders();
}

//name: runBoltz
//meta.cache: all
//meta.cache.invalidateOn: 0 0 1 * *
//input: string config
//output: string s
export async function runBoltz(config: string): Promise<string> {
  return await BoltzService.runBoltz(config);
}

//top-menu: Bio | Folding | Boltz-1...
//name: Folding
//input: dataframe df 
//input: column sequences {semType: Macromolecule}
//output: dataframe result
export async function folding(df: DG.DataFrame, sequences: DG.Column): Promise<DG.DataFrame> {
  return await BoltzService.folding(df, sequences);
}

//top-menu: Chem | Docking | Boltz-1...
//name: Docking
//input: dataframe df
//input: column molecules {semType: Molecule}
//input: string config {choices: Docking: getBoltzConfigFolders} [Folder with config files for docking]
//output: dataframe result
export async function docking(df: DG.DataFrame, molecules: DG.Column, config: string): Promise<DG.DataFrame> {
  return await BoltzService.docking(df, molecules, config);
}

//name: Boltz-1
//tags: panel, chem, widgets
//input: semantic_value molecule { semType: Molecule3D }
//condition: Boltz1:isApplicableBoltz(molecule)
//output: widget result
export async function boltzWidget(molecule: DG.SemanticValue): Promise<DG.Widget<any> | null> {
  return await BoltzService.boltzWidget(molecule);
}

//name: isApplicableBoltz
//input: string molecule
//output: bool result
export function isApplicableBoltz(molecule: string): boolean {
  return molecule.includes('confidence_score');
}

//tags: app
//name: Boltz-1 App
//output: view v
//meta.browsePath: Bio
export async function boltz1App(): Promise<DG.ViewBase> {
  return new Boltz1AppView().getView();
}