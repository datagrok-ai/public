import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Oligo Toolkit
//tags: app
//output: view result
//meta.icon: img/icons/toolkit.png
//meta.browsePath: Peptides | Oligo Toolkit
export async function oligoToolkitApp() {
  return PackageFunctions.oligoToolkitApp();
}

//tags: init
export async function init() {
  return PackageFunctions.init();
}

//name: Oligo Translator
//tags: app
//output: view result
//meta.icon: img/icons/translator.png
//meta.browsePath: Peptides | Oligo Toolkit
export async function oligoTranslatorApp() {
  return PackageFunctions.oligoTranslatorApp();
}

//name: Oligo Pattern
//tags: app
//output: view result
//meta.icon: img/icons/pattern.png
//meta.browsePath: Peptides | Oligo Toolkit
export async function oligoPatternApp() {
  return PackageFunctions.oligoPatternApp();
}

//name: Oligo Structure
//tags: app
//output: view result
//meta.icon: img/icons/structure.png
//meta.browsePath: Peptides | Oligo Toolkit
export async function oligoStructureApp() {
  return PackageFunctions.oligoStructureApp();
}

//name: getTranslationHelper
//output: object result
export async function getTranslationHelper() {
  return PackageFunctions.getTranslationHelper();
}

//name: getCodeToWeightsMap
//output: object result
export function getCodeToWeightsMap() {
  return PackageFunctions.getCodeToWeightsMap();
}

//name: validateSequence
//input: string sequence 
//output: bool result
export function validateSequence(sequence: string) {
  return PackageFunctions.validateSequence(sequence);
}

//name: validateSequence
//input: string sequence 
//input: bool invert 
//output: string result
export function getMolfileFromGcrsSequence(sequence: string, invert: boolean) {
  return PackageFunctions.getMolfileFromGcrsSequence(sequence, invert);
}

//name: linkStrands
//input: object strands 
//output: string result
export function linkStrands(strands: any) {
  return PackageFunctions.linkStrands(strands);
}

//name: demoOligoTranslator
//description: Translate oligonucleotide sequences across various formats accepted by different synthesizers
//meta.demoPath: Bioinformatics | Oligo Toolkit | Translator
//meta.path: /apps/Tutorials/Demo/Bioinformatics/Oligonucleotide%20Sequence:%20Translate
//meta.demoSkip: GROK-14320
export async function demoTranslateSequence() {
  return PackageFunctions.demoTranslateSequence();
}

//name: demoOligoPattern
//description: Design a modification pattern for an oligonucleotide sequence
//meta.demoPath: Bioinformatics | Oligo Toolkit | Pattern
//meta.path: %20/apps/Tutorials/Demo/Bioinformatics/Oligonucleotide%20Sequence:%20Visualize%20duplex
export async function demoOligoPattern() {
  return PackageFunctions.demoOligoPattern();
}

//name: demoOligoStructure
//description: Visualize duplex and save SDF
//meta.demoPath: Bioinformatics | Oligo Toolkit | Structure
//meta.path: %20/apps/Tutorials/Demo/Bioinformatics/Oligonucleotide%20Sequence:%20Visualize%20duplex
export async function demoOligoStructure() {
  return PackageFunctions.demoOligoStructure();
}

//name: translateOligonucleotideSequence
//input: string sequence 
//input: string sourceFormat 
//input: string targetFormat 
//output: string result
export async function translateOligonucleotideSequence(sequence: string, sourceFormat: string, targetFormat: string) {
  return PackageFunctions.translateOligonucleotideSequence(sequence, sourceFormat, targetFormat);
}

//name: polyToolConvert
//description: Perform cyclization of polymers
//top-menu: Bio | PolyTool | Convert...
export async function polyToolConvertTopMenu() {
  return PackageFunctions.polyToolConvertTopMenu();
}

//name: getPolyToolConvertEditor
//tags: editor
//input: funccall call 
//output: column result
export async function getPolyToolConvertEditor(call: DG.FuncCall) {
  return PackageFunctions.getPolyToolConvertEditor(call);
}

//name: polyToolConvert2
//input: dataframe table 
//input: column seqCol { caption: Sequence }
//input: bool generateHelm { default: true }
//input: bool chiralityEngine { default: true }
//input: object rules 
//output: column result
//editor: SequenceTranslator:getPolyToolConvertEditor
export async function polyToolConvert2(table: DG.DataFrame, seqCol: DG.Column, generateHelm: boolean, chiralityEngine: boolean, rules: string[]) {
  return PackageFunctions.polyToolConvert2(table, seqCol, generateHelm, chiralityEngine, rules);
}

//name: polyToolEnumerateHelm
//description: Perform cyclization of polymers
//top-menu: Bio | PolyTool | Enumerate HELM...
export async function polyToolEnumerateHelmTopMenu() {
  return PackageFunctions.polyToolEnumerateHelmTopMenu();
}

//name: polyToolEnumerateChem
//description: Perform cyclization of polymers
//top-menu: Bio | PolyTool | Enumerate Chem...
export async function polyToolEnumerateChemTopMenu() {
  return PackageFunctions.polyToolEnumerateChemTopMenu();
}

//name: polyToolColumnChoice
//input: dataframe df { description: Input data table }
//input: column macroMolecule 
export async function polyToolColumnChoice(df: DG.DataFrame, macroMolecule: DG.Column) {
  return PackageFunctions.polyToolColumnChoice(df, macroMolecule);
}

//name: createMonomerLibraryForPolyTool
//input: file file 
export async function createMonomerLibraryForPolyTool(file: DG.FileInfo) {
  return PackageFunctions.createMonomerLibraryForPolyTool(file);
}

//name: HELM Enumerator
//tags: app
//meta.icon: img/icons/structure.png
//meta.browsePath: Peptides | PolyTool
export async function ptEnumeratorHelmApp() {
  return PackageFunctions.ptEnumeratorHelmApp();
}

//name: Chem Enumerator
//tags: app
//meta.icon: img/icons/structure.png
//meta.browsePath: Peptides | PolyTool
export async function ptEnumeratorChemApp() {
  return PackageFunctions.ptEnumeratorChemApp();
}

//name: Polytool Helm Enumerator dialog
//input: object cell { nullable: true }
export async function getPtHelmEnumeratorDialog(cell?: any) {
  return PackageFunctions.getPtHelmEnumeratorDialog(cell);
}

//name: Polytool Chem Enumerator dialog
//input: object cell { nullable: true }
export async function getPtChemEnumeratorDialog(cell?: any) {
  return PackageFunctions.getPtChemEnumeratorDialog(cell);
}

//name: Combine Sequences
//top-menu: Bio | PolyTool | Combine Sequences...
export async function getPolyToolCombineDialog() {
  return PackageFunctions.getPolyToolCombineDialog();
}

//name: applyNotationProviderForHarmonizedSequence
//input: column col 
//input: string separator 
export function applyNotationProviderForCyclized(col: any, separator: string) {
  return PackageFunctions.applyNotationProviderForCyclized(col, separator);
}
