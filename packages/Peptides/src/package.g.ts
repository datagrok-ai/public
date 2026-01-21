import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//meta.role: init
export async function initPeptides() : Promise<void> {
  await PackageFunctions.initPeptides();
}

//output: view result
export function Peptides() : any {
  return PackageFunctions.Peptides();
}

//name: Bio Peptides
//top-menu: Bio | Analyze | SAR...
export function peptidesDialog() : any {
  return PackageFunctions.peptidesDialog();
}

//input: viewer v 
export async function testInitFunctionPeptides(v: any) : Promise<void> {
  await PackageFunctions.testInitFunctionPeptides(v);
}

//name: Peptides
//input: column col { semType: Macromolecule }
//output: widget result
//meta.role: widgets,panel
export function peptidesPanel(col: DG.Column) : any {
  return PackageFunctions.peptidesPanel(col);
}

//name: Sequence Variability Map
//description: Peptides Sequence Variability Map Viewer
//output: viewer result
//meta.icon: files/icons/peptide-sar-viewer.svg
//meta.role: viewer
export function monomerPosition() : any {
  return PackageFunctions.monomerPosition();
}

//name: Most Potent Residues
//description: Peptides Most Potent Residues Viewer
//output: viewer result
//meta.icon: files/icons/peptide-sar-vertical-viewer.svg
//meta.role: viewer
export function mostPotentResidues() : any {
  return PackageFunctions.mostPotentResidues();
}

//name: Sequence Mutation Cliffs
//description: Mutation Cliffs Line Chart
//output: viewer result
//meta.icon: files/icons/sequence-statistics-viewer.svg
//meta.role: viewer
export function mutationCliffs() : any {
  return PackageFunctions.mutationCliffs();
}

//name: Logo Summary Table
//output: viewer result
//meta.icon: files/icons/logo-summary-viewer.svg
//meta.role: viewer
export function logoSummaryTable() : any {
  return PackageFunctions.logoSummaryTable();
}

//name: Sequence Position Statistics
//output: viewer result
//meta.icon: files/icons/sequence-statistics-viewer.svg
//meta.role: viewer
export function sequencePositionStatistics() : any {
  return PackageFunctions.sequencePositionStatistics();
}

//name: Active peptide selection
//output: viewer result
//meta.role: viewer
export function clusterMaxActivity() : any {
  return PackageFunctions.clusterMaxActivity();
}

//name: Manual Alignment
//input: string _monomer { semType: Monomer }
//output: widget result
//meta.role: widgets,panel
export function manualAlignment(_monomer: string) : any {
  return PackageFunctions.manualAlignment(_monomer);
}

//name: Peptide SAR
//description: Peptide SAR Analysis demo on peptide sequences in FASTA format
//meta.demoPath: Bioinformatics | Peptide SAR
//meta.isDemoDashboard: true
export async function macromoleculeSarFastaDemo() : Promise<void> {
  await PackageFunctions.macromoleculeSarFastaDemo();
}

//name: LST Pie Chart
//output: grid_cell_renderer result
//meta.cellType: lst-pie-chart
//meta.gridChart: true
//meta.role: cellRenderer
export function lstPiechartCellRenderer() : any {
  return PackageFunctions.lstPiechartCellRenderer();
}
