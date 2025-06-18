import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace scripts {

}

export namespace funcs {
  export async function init(): Promise<any> {
    return await grok.functions.call('@datagrok/sequence-translator:Init', {});
  }

  export async function oligoToolkitApp(): Promise<any> {
    return await grok.functions.call('@datagrok/sequence-translator:OligoToolkitApp', {});
  }

  export async function oligoTranslatorApp(): Promise<any> {
    return await grok.functions.call('@datagrok/sequence-translator:OligoTranslatorApp', {});
  }

  export async function oligoPatternApp(): Promise<any> {
    return await grok.functions.call('@datagrok/sequence-translator:OligoPatternApp', {});
  }

  export async function oligoStructureApp(): Promise<any> {
    return await grok.functions.call('@datagrok/sequence-translator:OligoStructureApp', {});
  }

  export async function getTranslationHelper(): Promise<any> {
    return await grok.functions.call('@datagrok/sequence-translator:GetTranslationHelper', {});
  }

  export async function getCodeToWeightsMap(): Promise<any> {
    return await grok.functions.call('@datagrok/sequence-translator:GetCodeToWeightsMap', {});
  }

  export async function validateSequence(sequence: string): Promise<any> {
    return await grok.functions.call('@datagrok/sequence-translator:ValidateSequence', { sequence });
  }

  export async function getMolfileFromGcrsSequence(sequence: string, invert: boolean): Promise<any> {
    return await grok.functions.call('@datagrok/sequence-translator:GetMolfileFromGcrsSequence', { sequence, invert });
  }

  export async function linkStrands(strands: any): Promise<any> {
    return await grok.functions.call('@datagrok/sequence-translator:LinkStrands', { strands });
  }

  //Translate oligonucleotide sequences across various formats accepted by different synthesizers
  export async function demoTranslateSequence(): Promise<any> {
    return await grok.functions.call('@datagrok/sequence-translator:DemoTranslateSequence', {});
  }

  //Design a modification pattern for an oligonucleotide sequence
  export async function demoOligoPattern(): Promise<any> {
    return await grok.functions.call('@datagrok/sequence-translator:DemoOligoPattern', {});
  }

  //Visualize duplex and save SDF
  export async function demoOligoStructure(): Promise<any> {
    return await grok.functions.call('@datagrok/sequence-translator:DemoOligoStructure', {});
  }

  export async function translateOligonucleotideSequence(sequence: string, sourceFormat: string, targetFormat: string): Promise<any> {
    return await grok.functions.call('@datagrok/sequence-translator:TranslateOligonucleotideSequence', { sequence, sourceFormat, targetFormat });
  }

  //Perform cyclization of polymers
  export async function polyToolConvertTopMenu(): Promise<any> {
    return await grok.functions.call('@datagrok/sequence-translator:PolyToolConvertTopMenu', {});
  }

  export async function getPolyToolConvertEditor(call: any): Promise<any> {
    return await grok.functions.call('@datagrok/sequence-translator:GetPolyToolConvertEditor', { call });
  }

  export async function polyToolConvert2(table: DG.DataFrame, seqCol: DG.Column, generateHelm: boolean, chiralityEngine: boolean, rules: any): Promise<any> {
    return await grok.functions.call('@datagrok/sequence-translator:PolyToolConvert2', { table, seqCol, generateHelm, chiralityEngine, rules });
  }

  //Perform cyclization of polymers
  export async function polyToolEnumerateHelmTopMenu(): Promise<any> {
    return await grok.functions.call('@datagrok/sequence-translator:PolyToolEnumerateHelmTopMenu', {});
  }

  //Perform cyclization of polymers
  export async function polyToolEnumerateChemTopMenu(): Promise<any> {
    return await grok.functions.call('@datagrok/sequence-translator:PolyToolEnumerateChemTopMenu', {});
  }

  export async function polyToolColumnChoice(df: DG.DataFrame, macroMolecule: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/sequence-translator:PolyToolColumnChoice', { df, macroMolecule });
  }

  export async function createMonomerLibraryForPolyTool(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('@datagrok/sequence-translator:CreateMonomerLibraryForPolyTool', { file });
  }

  export async function ptEnumeratorHelmApp(): Promise<any> {
    return await grok.functions.call('@datagrok/sequence-translator:PtEnumeratorHelmApp', {});
  }

  export async function ptEnumeratorChemApp(): Promise<any> {
    return await grok.functions.call('@datagrok/sequence-translator:PtEnumeratorChemApp', {});
  }

  export async function getPtHelmEnumeratorDialog(cell: any): Promise<any> {
    return await grok.functions.call('@datagrok/sequence-translator:GetPtHelmEnumeratorDialog', { cell });
  }

  export async function getPtChemEnumeratorDialog(cell: any): Promise<any> {
    return await grok.functions.call('@datagrok/sequence-translator:GetPtChemEnumeratorDialog', { cell });
  }

  export async function getPolyToolCombineDialog(): Promise<any> {
    return await grok.functions.call('@datagrok/sequence-translator:GetPolyToolCombineDialog', {});
  }

  export async function applyNotationProviderForCyclized(col: DG.Column, separator: string): Promise<any> {
    return await grok.functions.call('@datagrok/sequence-translator:ApplyNotationProviderForCyclized', { col, separator });
  }
}
