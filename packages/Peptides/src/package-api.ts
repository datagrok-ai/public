import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace scripts {
  export async function smiTo3D(smiles: string): Promise<string> {
    return await grok.functions.call('@datagrok/peptides:SmiTo3D', { smiles });
  }
}

export namespace funcs {
  export async function initPeptides(): Promise<any> {
    return await grok.functions.call('@datagrok/peptides:InitPeptides', {});
  }

  export async function peptides(): Promise<any> {
    return await grok.functions.call('@datagrok/peptides:Peptides', {});
  }

  export async function peptidesDialog(): Promise<any> {
    return await grok.functions.call('@datagrok/peptides:PeptidesDialog', {});
  }

  export async function testInitFunctionPeptides(v: any): Promise<any> {
    return await grok.functions.call('@datagrok/peptides:TestInitFunctionPeptides', { v });
  }

  export async function peptidesPanel(col: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/peptides:PeptidesPanel', { col });
  }

  //Peptides Sequence Variability Map Viewer
  export async function monomerPosition(): Promise<any> {
    return await grok.functions.call('@datagrok/peptides:MonomerPosition', {});
  }

  //Peptides Most Potent Residues Viewer
  export async function mostPotentResidues(): Promise<any> {
    return await grok.functions.call('@datagrok/peptides:MostPotentResidues', {});
  }

  export async function logoSummaryTable(): Promise<any> {
    return await grok.functions.call('@datagrok/peptides:LogoSummaryTable', {});
  }

  export async function sequencePositionStatistics(): Promise<any> {
    return await grok.functions.call('@datagrok/peptides:SequencePositionStatistics', {});
  }

  export async function clusterMaxActivity(): Promise<any> {
    return await grok.functions.call('@datagrok/peptides:ClusterMaxActivity', {});
  }

  export async function manualAlignment(_monomer: string): Promise<any> {
    return await grok.functions.call('@datagrok/peptides:ManualAlignment', { _monomer });
  }

  //Peptide SAR Analysis demo on peptide sequences in FASTA format
  export async function macromoleculeSarFastaDemo(): Promise<any> {
    return await grok.functions.call('@datagrok/peptides:MacromoleculeSarFastaDemo', {});
  }

  export async function lstPiechartCellRenderer(): Promise<any> {
    return await grok.functions.call('@datagrok/peptides:LstPiechartCellRenderer', {});
  }
}
