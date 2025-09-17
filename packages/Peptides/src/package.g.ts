import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: initPeptides
//tags: init
export async function initPeptides() {
  return PackageFunctions.initPeptides();
}

//name: Peptides
//output: view result
export function Peptides() {
  return PackageFunctions.Peptides();
}

//name: Bio Peptides
//top-menu: Bio | Analyze | SAR...
export function peptidesDialog() {
  return PackageFunctions.peptidesDialog();
}

//name: testInitFunctionPeptides
//input: viewer v 
export async function testInitFunctionPeptides(v: any) {
  return PackageFunctions.testInitFunctionPeptides(v);
}

//name: Peptides
//tags: panel, widgets
//input: column col { semType: Macromolecule }
//output: widget result
export function peptidesPanel(col: DG.Column) {
  return PackageFunctions.peptidesPanel(col);
}

//name: Sequence Variability Map
//description: Peptides Sequence Variability Map Viewer
//tags: viewer
//output: viewer result
//meta.icon: files/icons/peptide-sar-viewer.svg
export function monomerPosition() {
  return PackageFunctions.monomerPosition();
}

//name: Most Potent Residues
//description: Peptides Most Potent Residues Viewer
//tags: viewer
//output: viewer result
//meta.icon: files/icons/peptide-sar-vertical-viewer.svg
export function mostPotentResidues() {
  return PackageFunctions.mostPotentResidues();
}

//name: Logo Summary Table
//tags: viewer
//output: viewer result
//meta.icon: files/icons/logo-summary-viewer.svg
export function logoSummaryTable() {
  return PackageFunctions.logoSummaryTable();
}

//name: Sequence Position Statistics
//tags: viewer
//output: viewer result
//meta.icon: files/icons/sequence-statistics-viewer.svg
export function sequencePositionStatistics() {
  return PackageFunctions.sequencePositionStatistics();
}

//name: Active peptide selection
//tags: viewer
//output: viewer result
export function clusterMaxActivity() {
  return PackageFunctions.clusterMaxActivity();
}

//name: Manual Alignment
//tags: panel, widgets
//input: string _monomer { semType: Monomer }
//output: widget result
export function manualAlignment(_monomer: string) {
  return PackageFunctions.manualAlignment(_monomer);
}

//name: Peptide SAR
//description: Peptide SAR Analysis demo on peptide sequences in FASTA format
//meta.demoPath: Bioinformatics | Peptide SAR
//meta.demoSkip: GROK-14320
export async function macromoleculeSarFastaDemo() {
  return PackageFunctions.macromoleculeSarFastaDemo();
}

//name: LST Pie Chart
//tags: cellRenderer
//output: grid_cell_renderer result
//meta.cellType: lst-pie-chart
//meta.gridChart: true
export function lstPiechartCellRenderer() {
  return PackageFunctions.lstPiechartCellRenderer();
}
