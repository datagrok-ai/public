import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//meta.role: init
export async function init() : Promise<void> {
  await PackageFunctions.init();
}

//meta.role: autostart
//meta.autostartImmediate: true
export async function autostart() : Promise<void> {
  await PackageFunctions.autostart();
}

//name: Demo App
//meta.role: adminApp,app
//meta.browsePath: Chem | Demo
//meta.icon: files/icons/demo.svg
//output: view result
export function demoApp() : any {
  return PackageFunctions.demoApp();
}

//name: Demo Viewer
//description: A viewer that is also a panel
//meta.role: viewer,panel
//meta.trellisable: true
//meta.icon: files/icons/viewer.svg
//output: viewer result
export function demoViewer() : any {
  return PackageFunctions.demoViewer();
}

//name: Molecule Panel
//meta.role: widgets,panel
//condition: true
//input: string smiles { semType: Molecule }
//output: widget result
export function moleculePanel(smiles: string) : any {
  return PackageFunctions.moleculePanel(smiles);
}

//name: Molecule Renderer
//meta.role: cellRenderer
//meta.cellType: Molecule
//meta.columnTags: quality=Molecule, foo=bar
//output: grid_cell_renderer result
export function moleculeRenderer() : any {
  return PackageFunctions.moleculeRenderer();
}

//name: Substructure Filter
//meta.role: filter
//meta.semType: Molecule
//meta.primaryFilter: true
//output: filter result
export function substructureFilter() : any {
  return PackageFunctions.substructureFilter();
}

//name: Import SDF
//meta.role: fileHandler
//meta.ext: sdf,mol
//input: string bytes
//output: list tables
export function importSdf(bytes: string) : any {
  return PackageFunctions.importSdf(bytes);
}

//name: Preview MOL
//meta.role: fileViewer
//meta.fileViewer: mol,mol2
//meta.fileViewerCheck: Demo:checkMol
//input: file file
//output: view result
export function previewMol(file: any) : any {
  return PackageFunctions.previewMol(file);
}

//name: Broken Renderer
//meta.role: cellRenderer
//output: grid_cell_renderer result
export function brokenRenderer() : any {
  return PackageFunctions.brokenRenderer();
}

//name: Column Editor
//meta.role: editor
//editor-for: AddNewColumn
//input: funccall call
export function columnEditor(call: any) : void {
  PackageFunctions.columnEditor(call);
}

//name: Demo Handler
//meta.role: scriptHandler
//meta.scriptHandler.language: demo
//meta.scriptHandler.extensions: dm
//meta.scriptHandler.commentStart: #
export function demoHandler() : void {
  PackageFunctions.demoHandler();
}

//name: To HELM
//description: Converts a sequence to HELM
//feature: domains/bio
//top-menu: Bio | Convert | To HELM...
//help-url: https://datagrok.ai/help/domains/bio
//meta.cache: true
//meta.cache.invalidateOn: 0 0 * * *
//meta.demoPath: Bioinformatics | To HELM
//meta.vectorFunc: true
//input: dataframe table
//input: column sequence { semType: Macromolecule; units: fasta }
//input: string method { choices: ["fast","slow"] } [The conversion method]
//input: int n = 5 { optional: true }
//output: column result { semType: Macromolecule }
export async function toHelm(table: DG.DataFrame, sequence: DG.Column, method: string, n: number) : Promise<any> {
  return await PackageFunctions.toHelm(table, sequence, method, n);
}

//input: string mol
//output: string molfile { semType: Molecule }
export function toMolfile(mol: string) : string {
  return PackageFunctions.toMolfile(mol);
}
