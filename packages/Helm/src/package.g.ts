import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//tags: init
//meta.role: init
export async function initHelm() : Promise<void> {
  await PackageFunctions.initHelm();
}

//tags: cellRenderer
//output: grid_cell_renderer result
//meta.columnTags: quality=Macromolecule, units=helm
//meta.cellType: helm
//meta.role: cellRenderer
export function helmCellRenderer() : any {
  return PackageFunctions.helmCellRenderer();
}

//description: Macromolecule
//tags: cellEditor
//input: grid_cell cell 
//meta.columnTags: quality=Macromolecule, units=helm
//meta.role: cellEditor
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
//tags: panel, widgets, bio
//input: semantic_value sequence { semType: Macromolecule }
//output: widget result
//meta.role: widgets,panel
//meta.domain: bio
export async function propertiesWidget(sequence: DG.SemanticValue) : Promise<any> {
  return await PackageFunctions.propertiesWidget(sequence);
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
//meta.role: valueEditor
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
