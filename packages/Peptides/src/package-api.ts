import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Scripts {
  export async function smiTo3D(smiles: string): Promise<string> {
    return await grok.functions.call('Peptides:SmiTo3D', { smiles });
  }
}

export namespace Funcs {
  export async function initPeptides(): Promise<any> {
    return await grok.functions.call('Peptides:InitPeptides', {});
  }

  export async function peptides(): Promise<any> {
    return await grok.functions.call('Peptides:Peptides', {});
  }

  export async function peptidesDialog(): Promise<any> {
    return await grok.functions.call('Peptides:PeptidesDialog', {});
  }

  export async function testInitFunctionPeptides(v: any): Promise<any> {
    return await grok.functions.call('Peptides:TestInitFunctionPeptides', { v });
  }

  export async function peptidesPanel(col: DG.Column): Promise<any> {
    return await grok.functions.call('Peptides:PeptidesPanel', { col });
  }

  //Peptides Sequence Variability Map Viewer
  export async function monomerPosition(): Promise<any> {
    return await grok.functions.call('Peptides:MonomerPosition', {});
  }

  //Peptides Most Potent Residues Viewer
  export async function mostPotentResidues(): Promise<any> {
    return await grok.functions.call('Peptides:MostPotentResidues', {});
  }

  export async function logoSummaryTable(): Promise<any> {
    return await grok.functions.call('Peptides:LogoSummaryTable', {});
  }

  export async function sequencePositionStatistics(): Promise<any> {
    return await grok.functions.call('Peptides:SequencePositionStatistics', {});
  }

  export async function clusterMaxActivity(): Promise<any> {
    return await grok.functions.call('Peptides:ClusterMaxActivity', {});
  }

  export async function manualAlignment(_monomer: string): Promise<any> {
    return await grok.functions.call('Peptides:ManualAlignment', { _monomer });
  }

  //Peptide SAR Analysis demo on peptide sequences in FASTA format
  export async function macromoleculeSarFastaDemo(): Promise<any> {
    return await grok.functions.call('Peptides:MacromoleculeSarFastaDemo', {});
  }

  export async function lstPiechartCellRenderer(): Promise<any> {
    return await grok.functions.call('Peptides:LstPiechartCellRenderer', {});
  }
}
