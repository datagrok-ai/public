import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//tags: init
export async function initHelm() : Promise<void> {
  await PackageFunctions.initHelm();
}

//description: Helm renderer service
//output: object result
export async function getHelmService() : Promise<any> {
  return await PackageFunctions.getHelmService();
}

//tags: cellRenderer
//output: grid_cell_renderer result
//meta.columnTags: quality=Macromolecule, units=helm
//meta.cellType: helm
export function helmCellRenderer() : any {
  return PackageFunctions.helmCellRenderer();
}

//description: Macromolecule
//tags: cellEditor
//input: grid_cell cell 
//meta.columnTags: quality=Macromolecule, units=helm
export function editMoleculeCell(cell: any) : void {
  PackageFunctions.editMoleculeCell(cell);
}

//name: Edit Helm...
//description: Adds editor
//input: semantic_value mol { semType: Macromolecule }
//meta.action: Edit Helm...
export function openEditor(mol: DG.SemanticValue) : void {
  PackageFunctions.openEditor(mol);
}

//name: Properties
//tags: panel, bio, helm, widgets
//input: semantic_value sequence { semType: Macromolecule }
//output: widget result
export function propertiesWidget(sequence: DG.SemanticValue) : any {
  return PackageFunctions.propertiesWidget(sequence);
}

//input: column col { semType: Macromolecule }
//output: column result
export function getMolfiles(col: DG.Column<any>) : any {
  return PackageFunctions.getMolfiles(col);
}

//tags: valueEditor
//input: string name { optional: true }
//input: object options { optional: true }
//output: object result
//meta.propertyType: string
//meta.semType: Macromolecule
export function helmInput(name: string, options: any) : any {
  return PackageFunctions.helmInput(name, options);
}

//output: object result
export function getHelmHelper() : any {
  return PackageFunctions.getHelmHelper();
}

//name: measureCellRenderer
export async function measureCellRenderer() : Promise<void> {
  await PackageFunctions.measureCellRenderer();
}

//name: Highlight Monomers
export async function highlightMonomers() : Promise<void> {
  await PackageFunctions.highlightMonomers();
}
