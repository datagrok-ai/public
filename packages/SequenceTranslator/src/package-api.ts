import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace scripts {

}

export namespace funcs {
  export async function init(): Promise<any> {
    return await grok.functions.call('SequenceTranslator:Init', {});
  }

  export async function oligoToolkitApp(): Promise<any> {
    return await grok.functions.call('SequenceTranslator:OligoToolkitApp', {});
  }

  export async function oligoTranslatorApp(): Promise<any> {
    return await grok.functions.call('SequenceTranslator:OligoTranslatorApp', {});
  }

  export async function oligoPatternApp(): Promise<any> {
    return await grok.functions.call('SequenceTranslator:OligoPatternApp', {});
  }

  export async function oligoStructureApp(): Promise<any> {
    return await grok.functions.call('SequenceTranslator:OligoStructureApp', {});
  }

  export async function getTranslationHelper(): Promise<any> {
    return await grok.functions.call('SequenceTranslator:GetTranslationHelper', {});
  }

  export async function getCodeToWeightsMap(): Promise<any> {
    return await grok.functions.call('SequenceTranslator:GetCodeToWeightsMap', {});
  }

  export async function validateSequence(sequence: string): Promise<any> {
    return await grok.functions.call('SequenceTranslator:ValidateSequence', { sequence });
  }

  export async function getMolfileFromGcrsSequence(sequence: string, invert: boolean): Promise<any> {
    return await grok.functions.call('SequenceTranslator:GetMolfileFromGcrsSequence', { sequence, invert });
  }

  export async function linkStrands(strands: any): Promise<any> {
    return await grok.functions.call('SequenceTranslator:LinkStrands', { strands });
  }

  //Translate oligonucleotide sequences across various formats accepted by different synthesizers
  export async function demoTranslateSequence(): Promise<any> {
    return await grok.functions.call('SequenceTranslator:DemoTranslateSequence', {});
  }

  //Design a modification pattern for an oligonucleotide sequence
  export async function demoOligoPattern(): Promise<any> {
    return await grok.functions.call('SequenceTranslator:DemoOligoPattern', {});
  }

  //Visualize duplex and save SDF
  export async function demoOligoStructure(): Promise<any> {
    return await grok.functions.call('SequenceTranslator:DemoOligoStructure', {});
  }

  export async function translateOligonucleotideSequence(sequence: string, sourceFormat: string, targetFormat: string): Promise<any> {
    return await grok.functions.call('SequenceTranslator:TranslateOligonucleotideSequence', { sequence, sourceFormat, targetFormat });
  }

  //Perform cyclization of polymers
  export async function polyToolConvertTopMenu(): Promise<any> {
    return await grok.functions.call('SequenceTranslator:PolyToolConvertTopMenu', {});
  }

  export async function getPolyToolConvertEditor(call: any): Promise<any> {
    return await grok.functions.call('SequenceTranslator:GetPolyToolConvertEditor', { call });
  }

  export async function polyToolConvert2(table: DG.DataFrame, seqCol: DG.Column, generateHelm: boolean, chiralityEngine: boolean, rules: any): Promise<any> {
    return await grok.functions.call('SequenceTranslator:PolyToolConvert2', { table, seqCol, generateHelm, chiralityEngine, rules });
  }

  //Perform cyclization of polymers
  export async function polyToolEnumerateHelmTopMenu(): Promise<any> {
    return await grok.functions.call('SequenceTranslator:PolyToolEnumerateHelmTopMenu', {});
  }

  //Perform cyclization of polymers
  export async function polyToolEnumerateChemTopMenu(): Promise<any> {
    return await grok.functions.call('SequenceTranslator:PolyToolEnumerateChemTopMenu', {});
  }

  export async function polyToolColumnChoice(df: DG.DataFrame, macroMolecule: DG.Column): Promise<any> {
    return await grok.functions.call('SequenceTranslator:PolyToolColumnChoice', { df, macroMolecule });
  }

  export async function createMonomerLibraryForPolyTool(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('SequenceTranslator:CreateMonomerLibraryForPolyTool', { file });
  }

  export async function ptEnumeratorHelmApp(): Promise<any> {
    return await grok.functions.call('SequenceTranslator:PtEnumeratorHelmApp', {});
  }

  export async function ptEnumeratorChemApp(): Promise<any> {
    return await grok.functions.call('SequenceTranslator:PtEnumeratorChemApp', {});
  }

  export async function getPtHelmEnumeratorDialog(cell: any): Promise<any> {
    return await grok.functions.call('SequenceTranslator:GetPtHelmEnumeratorDialog', { cell });
  }

  export async function getPtChemEnumeratorDialog(cell: any): Promise<any> {
    return await grok.functions.call('SequenceTranslator:GetPtChemEnumeratorDialog', { cell });
  }

  export async function getPolyToolCombineDialog(): Promise<any> {
    return await grok.functions.call('SequenceTranslator:GetPolyToolCombineDialog', {});
  }

  export async function applyNotationProviderForCyclized(col: DG.Column, separator: string): Promise<any> {
    return await grok.functions.call('SequenceTranslator:ApplyNotationProviderForCyclized', { col, separator });
  }
}
