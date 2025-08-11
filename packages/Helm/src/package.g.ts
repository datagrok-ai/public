import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: initHelm
//tags: init
export async function initHelm() {
  return PackageFunctions.initHelm();
}

//name: getHelmService
//description: Helm renderer service
//output: object result
export async function getHelmService() {
  return PackageFunctions.getHelmService();
}

//name: helmCellRenderer
//tags: cellRenderer
//output: grid_cell_renderer result
//meta.columnTags: quality=Macromolecule, units=helm
//meta.cellType: helm
export function helmCellRenderer() {
  return PackageFunctions.helmCellRenderer();
}

//name: editMoleculeCell
//description: Macromolecule
//tags: cellEditor
//input: grid_cell cell 
//meta.columnTags: quality=Macromolecule, units=helm
export function editMoleculeCell(cell: any) {
  return PackageFunctions.editMoleculeCell(cell);
}

//name: Edit Helm...
//description: Adds editor
//input: semantic_value mol { semType: Macromolecule }
//meta.action: Edit Helm...
export function openEditor(mol: DG.SemanticValue) {
  return PackageFunctions.openEditor(mol);
}

//name: Properties
//tags: panel, bio, helm, widgets
//input: semantic_value sequence { semType: Macromolecule }
//output: widget result
export function propertiesWidget(sequence: DG.SemanticValue) {
  return PackageFunctions.propertiesWidget(sequence);
}

//name: getMolfiles
//input: column col { semType: Macromolecule }
//output: column result
export function getMolfiles(col: any) {
  return PackageFunctions.getMolfiles(col);
}

//name: helmInput
//tags: valueEditor
//input: string name { optional: true; default: undefined }
//input: object options { optional: true; default: undefined }
//output: object result
//meta.propertyType: string
//meta.semType: Macromolecule
export function helmInput(name: string, options: any) {
  return PackageFunctions.helmInput(name, options);
}

//name: getHelmHelper
//output: object result
export function getHelmHelper() {
  return PackageFunctions.getHelmHelper();
}

//name: measureCellRenderer
export async function measureCellRenderer() {
  return PackageFunctions.measureCellRenderer();
}

//name: Highlight Monomers
export async function highlightMonomers() {
  return PackageFunctions.highlightMonomers();
}
