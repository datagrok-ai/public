import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Oligo Toolkit
//tags: app
//output: view result
//meta.icon: img/icons/toolkit.png
//meta.browsePath: Peptides | Oligo Toolkit
export async function oligoToolkitApp() : Promise<any> {
  return await PackageFunctions.oligoToolkitApp();
}

//tags: init
export async function init() : Promise<void> {
  await PackageFunctions.init();
}

//name: Oligo Translator
//tags: app
//output: view result
//meta.icon: img/icons/translator.png
//meta.browsePath: Peptides | Oligo Toolkit
export async function oligoTranslatorApp() : Promise<any> {
  return await PackageFunctions.oligoTranslatorApp();
}

//name: Oligo Pattern
//tags: app
//output: view result
//meta.icon: img/icons/pattern.png
//meta.browsePath: Peptides | Oligo Toolkit
export async function oligoPatternApp() : Promise<any> {
  return await PackageFunctions.oligoPatternApp();
}

//name: Oligo Structure
//tags: app
//output: view result
//meta.icon: img/icons/structure.png
//meta.browsePath: Peptides | Oligo Toolkit
export async function oligoStructureApp() : Promise<any> {
  return await PackageFunctions.oligoStructureApp();
}

//output: object result
export async function getTranslationHelper() : Promise<any> {
  return await PackageFunctions.getTranslationHelper();
}

//output: object result
export function getCodeToWeightsMap() : any {
  return PackageFunctions.getCodeToWeightsMap();
}

//input: string sequence 
//output: bool result
export function validateSequence(sequence: string) : boolean {
  return PackageFunctions.validateSequence(sequence);
}

//name: validateSequence
//input: string sequence 
//input: bool invert 
//output: string result
export function getMolfileFromGcrsSequence(sequence: string, invert: boolean) : string {
  return PackageFunctions.getMolfileFromGcrsSequence(sequence, invert);
}

//input: object strands 
//output: string result
export function linkStrands(strands: any) : string {
  return PackageFunctions.linkStrands(strands);
}

//name: demoOligoTranslator
//description: Translate oligonucleotide sequences across various formats accepted by different synthesizers
//meta.demoPath: Bioinformatics | Oligo Toolkit | Translator
//meta.path: /apps/Tutorials/Demo/Bioinformatics/Oligonucleotide%20Sequence:%20Translate
//meta.demoSkip: GROK-14320
export async function demoTranslateSequence() : Promise<void> {
  await PackageFunctions.demoTranslateSequence();
}

//description: Design a modification pattern for an oligonucleotide sequence
//meta.demoPath: Bioinformatics | Oligo Toolkit | Pattern
//meta.path: %20/apps/Tutorials/Demo/Bioinformatics/Oligonucleotide%20Sequence:%20Visualize%20duplex
export async function demoOligoPattern() : Promise<void> {
  await PackageFunctions.demoOligoPattern();
}

//description: Visualize duplex and save SDF
//meta.demoPath: Bioinformatics | Oligo Toolkit | Structure
//meta.path: %20/apps/Tutorials/Demo/Bioinformatics/Oligonucleotide%20Sequence:%20Visualize%20duplex
export async function demoOligoStructure() : Promise<void> {
  await PackageFunctions.demoOligoStructure();
}

//input: string sequence 
//input: string sourceFormat 
//input: string targetFormat 
//output: string result
export async function translateOligonucleotideSequence(sequence: string, sourceFormat: string, targetFormat: string) : Promise<string> {
  return await PackageFunctions.translateOligonucleotideSequence(sequence, sourceFormat, targetFormat);
}

//name: polyToolConvert
//description: Perform cyclization of polymers
//top-menu: Bio | PolyTool | Convert...
export async function polyToolConvertTopMenu() : Promise<void> {
  await PackageFunctions.polyToolConvertTopMenu();
}

//tags: editor
//input: funccall call 
//output: column result
export async function getPolyToolConvertEditor(call: DG.FuncCall) : Promise<any> {
  return await PackageFunctions.getPolyToolConvertEditor(call);
}

//input: dataframe table 
//input: column seqCol { caption: Sequence }
//input: bool generateHelm { default: true }
//input: bool chiralityEngine { default: true }
//input: object rules 
//output: column result
//editor: SequenceTranslator:getPolyToolConvertEditor
export async function polyToolConvert2(table: DG.DataFrame, seqCol: DG.Column, generateHelm: boolean, chiralityEngine: boolean, rules: string[]) : Promise<any> {
  return await PackageFunctions.polyToolConvert2(table, seqCol, generateHelm, chiralityEngine, rules);
}

//name: polyToolEnumerateHelm
//description: Perform cyclization of polymers
//top-menu: Bio | PolyTool | Enumerate HELM...
export async function polyToolEnumerateHelmTopMenu() : Promise<void> {
  await PackageFunctions.polyToolEnumerateHelmTopMenu();
}

//name: polyToolEnumerateChem
//description: Perform cyclization of polymers
//top-menu: Bio | PolyTool | Enumerate Chem...
export async function polyToolEnumerateChemTopMenu() : Promise<void> {
  await PackageFunctions.polyToolEnumerateChemTopMenu();
}

//input: dataframe df { description: Input data table }
//input: column macroMolecule 
export async function polyToolColumnChoice(df: DG.DataFrame, macroMolecule: DG.Column) : Promise<void> {
  await PackageFunctions.polyToolColumnChoice(df, macroMolecule);
}

//input: file file 
export async function createMonomerLibraryForPolyTool(file: DG.FileInfo) : Promise<void> {
  await PackageFunctions.createMonomerLibraryForPolyTool(file);
}

//name: HELM Enumerator
//tags: app
//meta.icon: img/icons/structure.png
//meta.browsePath: Peptides | PolyTool
export async function ptEnumeratorHelmApp() : Promise<void> {
  await PackageFunctions.ptEnumeratorHelmApp();
}

//name: Chem Enumerator
//tags: app
//meta.icon: img/icons/structure.png
//meta.browsePath: Peptides | PolyTool
export async function ptEnumeratorChemApp() : Promise<void> {
  await PackageFunctions.ptEnumeratorChemApp();
}

//name: Polytool Helm Enumerator dialog
//input: object cell { nullable: true }
export async function getPtHelmEnumeratorDialog(cell?: any) : Promise<void> {
  await PackageFunctions.getPtHelmEnumeratorDialog(cell);
}

//name: Polytool Chem Enumerator dialog
//input: object cell { nullable: true }
export async function getPtChemEnumeratorDialog(cell?: any) : Promise<void> {
  await PackageFunctions.getPtChemEnumeratorDialog(cell);
}

//name: Combine Sequences
//top-menu: Bio | PolyTool | Combine Sequences...
export async function getPolyToolCombineDialog() : Promise<void> {
  await PackageFunctions.getPolyToolCombineDialog();
}

//name: applyNotationProviderForHarmonizedSequence
//input: column col 
//input: string separator 
export function applyNotationProviderForCyclized(col: DG.Column<any>, separator: string) : void {
  PackageFunctions.applyNotationProviderForCyclized(col, separator);
}
