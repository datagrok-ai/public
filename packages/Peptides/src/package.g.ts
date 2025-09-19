import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//tags: init
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
//tags: panel, widgets
//input: column col { semType: Macromolecule }
//output: widget result
export function peptidesPanel(col: DG.Column) : any {
  return PackageFunctions.peptidesPanel(col);
}

//name: Sequence Variability Map
//description: Peptides Sequence Variability Map Viewer
//tags: viewer
//output: viewer result
//meta.icon: files/icons/peptide-sar-viewer.svg
export function monomerPosition() : any {
  return PackageFunctions.monomerPosition();
}

//name: Most Potent Residues
//description: Peptides Most Potent Residues Viewer
//tags: viewer
//output: viewer result
//meta.icon: files/icons/peptide-sar-vertical-viewer.svg
export function mostPotentResidues() : any {
  return PackageFunctions.mostPotentResidues();
}

//name: Logo Summary Table
//tags: viewer
//output: viewer result
//meta.icon: files/icons/logo-summary-viewer.svg
export function logoSummaryTable() : any {
  return PackageFunctions.logoSummaryTable();
}

//name: Sequence Position Statistics
//tags: viewer
//output: viewer result
//meta.icon: files/icons/sequence-statistics-viewer.svg
export function sequencePositionStatistics() : any {
  return PackageFunctions.sequencePositionStatistics();
}

//name: Active peptide selection
//tags: viewer
//output: viewer result
export function clusterMaxActivity() : any {
  return PackageFunctions.clusterMaxActivity();
}

//name: Manual Alignment
//tags: panel, widgets
//input: string _monomer { semType: Monomer }
//output: widget result
export function manualAlignment(_monomer: string) : any {
  return PackageFunctions.manualAlignment(_monomer);
}

//name: Peptide SAR
//description: Peptide SAR Analysis demo on peptide sequences in FASTA format
//meta.demoPath: Bioinformatics | Peptide SAR
export async function macromoleculeSarFastaDemo() : Promise<void> {
  await PackageFunctions.macromoleculeSarFastaDemo();
}

//name: LST Pie Chart
//tags: cellRenderer
//output: grid_cell_renderer result
//meta.cellType: lst-pie-chart
//meta.gridChart: true
export function lstPiechartCellRenderer() : any {
  return PackageFunctions.lstPiechartCellRenderer();
}
