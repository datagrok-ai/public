import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule/consts';
import {SeqTemps} from '@datagrok-libraries/bio/src/utils/macromolecule/seq-handler';

import {OligoToolkitPackage} from './apps/common/model/oligo-toolkit-package';
import {APP_NAME} from './apps/common/view/const';
import {getSpecifiedAppUI} from './apps/common/view/utils';
import {CombinedAppUI} from './apps/common/view/combined-app-ui';
import {linkStrandsV3000} from './apps/structure/model/mol-transformations';
import {SequenceToMolfileConverter} from './apps/structure/model/sequence-to-molfile';
import {demoOligoPatternUI, demoOligoStructureUI, demoOligoTranslatorUI} from './demo/demo-st-ui';
import {getExternalAppViewFactories} from './plugins/mermade';
import {defaultErrorHandler} from './utils/err-info';

//polytool specific
import {polyToolConvert, polyToolConvertUI} from './polytool/pt-dialog';
import {polyToolEnumerateChemUI} from './polytool/pt-dialog';
import {polyToolEnumerateHelmUI} from './polytool/pt-enumerate-seq-dialog';
import {_setPeptideColumn} from './polytool/utils';
import {PolyToolCsvLibHandler} from './polytool/csv-to-json-monomer-lib-converter';
import {ITranslationHelper} from './types';
import {PolyToolConvertFuncEditor} from './polytool/pt-convert-editor';
import {CyclizedNotationProvider} from './utils/cyclized';
import {getSeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {PolyToolTags} from './consts';
import {getHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import { getPTCombineDialog } from './polytool/pt-combine-dialog';

export const _package: OligoToolkitPackage = new OligoToolkitPackage({debug: true}/**/);

let initSequenceTranslatorPromise: Promise<void> | null = null;

//tags: init
export async function init(): Promise<void> {
  if (initSequenceTranslatorPromise === null)
    _package.startInit(initSequenceTranslatorPromise = initSequenceTranslatorInt());

  return initSequenceTranslatorPromise;
}

async function initSequenceTranslatorInt(): Promise<void> {
  const [helmHelper] = await Promise.all([
    getHelmHelper(),
  ]);
  _package.completeInit(helmHelper);
}

//name: Oligo Toolkit
//meta.icon: img/icons/toolkit.png
//meta.browsePath: Peptides | Oligo Toolkit
//tags: app
//output: view v
export async function oligoToolkitApp(): Promise<DG.ViewBase> {
  await _package.initLibData();
  const externalViewFactories = await getExternalAppViewFactories(_package);
  if (!externalViewFactories)
    throw new Error('External app view factories not loaded');
  const appUI = new CombinedAppUI(externalViewFactories!, _package);
  const view = await appUI.getAppView();
  return view;
}

//name: Oligo Translator
//meta.icon: img/icons/translator.png
//meta.browsePath: Peptides | Oligo Toolkit
//tags: app
//output: view v
export async function oligoTranslatorApp(): Promise<DG.ViewBase> {
  const view = await getSpecifiedAppView(APP_NAME.TRANSLATOR);
  return view;
}

//name: Oligo Pattern
//meta.icon: img/icons/pattern.png
//meta.browsePath: Peptides | Oligo Toolkit
//tags: app
//output: view v
export async function oligoPatternApp(): Promise<DG.ViewBase> {
  const view = await getSpecifiedAppView(APP_NAME.PATTERN);
  return view;
}

//name: Oligo Structure
//meta.icon: img/icons/structure.png
//meta.browsePath: Peptides | Oligo Toolkit
//tags: app
//output: view v
export async function oligoStructureApp(): Promise<DG.ViewBase> {
  const view = await getSpecifiedAppView(APP_NAME.STRUCTURE);
  return view;
}

//name: getTranslationHelper
//output: object result
export async function getTranslationHelper(): Promise<ITranslationHelper> {
  await _package.initLibData();
  return _package;
}

//name: getCodeToWeightsMap
//output: object result
export function getCodeToWeightsMap(): { [key: string]: number } {
  const monomerLibWrapper = _package.monomerLibWrapper;
  const map = monomerLibWrapper.getCodesToWeightsMap();
  return Object.fromEntries(map);
}

//name: validateSequence
//input: string sequence
//output: bool result
export function validateSequence(sequence: string): boolean {
  const validator = _package.createSequenceValidator(sequence);
  const format = _package.createFormatDetector(sequence).getFormat();
  return (format === null) ? false : validator.isValidSequence(format!);
}

//name: validateSequence
//input: string sequence
//input: bool invert
//output: string result
export function getMolfileFromGcrsSequence(sequence: string, invert: boolean): string {
  return (new SequenceToMolfileConverter(sequence, invert, 'GCRS')).convert();
}

//name: linkStrands
//input: object strands
//output: string result
export function linkStrands(strands: { senseStrands: string[], antiStrands: string[] }): string {
  return linkStrandsV3000(strands, true);
}

//name: demoOligoTranslator
//meta.demoPath: Bioinformatics | Oligo Toolkit | Translator
//description: Translate oligonucleotide sequences across various formats accepted by different synthesizers
//meta.path: /apps/Tutorials/Demo/Bioinformatics/Oligonucleotide%20Sequence:%20Translate
//meta.demoSkip: GROK-14320
export async function demoTranslateSequence(): Promise<void> {
  await demoOligoTranslatorUI();
}

//name: demoOligoPattern
//meta.demoPath: Bioinformatics | Oligo Toolkit | Pattern
//description: Design a modification pattern for an oligonucleotide sequence
//meta.path:%20/apps/Tutorials/Demo/Bioinformatics/Oligonucleotide%20Sequence:%20Visualize%20duplex
export async function demoOligoPattern(): Promise<void> {
  await demoOligoPatternUI();
}

//name: demoOligoStructure
//meta.demoPath: Bioinformatics | Oligo Toolkit | Structure
//description: Visualize duplex and save SDF
//meta.path:%20/apps/Tutorials/Demo/Bioinformatics/Oligonucleotide%20Sequence:%20Visualize%20duplex
export async function demoOligoStructure(): Promise<void> {
  await demoOligoStructureUI();
}

//name: translateOligonucleotideSequence
//input: string sequence
//input: string sourceFormat
//input: string targetFormat
//output: string result
export async function translateOligonucleotideSequence(
  sequence: string, sourceFormat: string, targetFormat: string
): Promise<string> {
  await _package.initLibData();
  return _package.createFormatConverter(sequence, sourceFormat).convertTo(targetFormat);
}

async function getSpecifiedAppView(appName: string): Promise<DG.ViewBase> {
  await _package.initLibData();
  const appUI = getSpecifiedAppUI(appName, _package);
  const view = await appUI.getAppView();
  return view;
}

//top-menu: Bio | PolyTool | Convert...
//name: polyToolConvert
//description: Perform cyclization of polymers
export async function polyToolConvertTopMenu(): Promise<void> {
  await polyToolConvertUI();
}

// //top-menu: Bio | PolyTool | Unrule...
// //name: polyToolUnrule
// //description: Perform uncyclization of polymers by rules
// export async function polyToolUnruleTopMenu(): Promise<void> {
//   await polyToolUnruleUI();
// }

//name: getPolyToolConvertEditor
//tags: editor
//input: funccall call
//output: column resCol
export async function getPolyToolConvertEditor(call: DG.FuncCall): Promise<DG.Column<string> | null> {
  const funcEditor = await PolyToolConvertFuncEditor.create(call);
  return await funcEditor.showDialog();
}

//name: polyToolConvert2
//input: dataframe table
//input: column seqCol { caption: Sequence }
//input: bool generateHelm = true
//input: bool chiralityEngine = true
//input: object rules
//output: column resCol
//editor: SequenceTranslator:getPolyToolConvertEditor
export async function polyToolConvert2(table: DG.DataFrame,
  seqCol: DG.Column, generateHelm: boolean, chiralityEngine: boolean, rules: string[]
): Promise<DG.Column<string>> {
  const ptConvertRes = await polyToolConvert(seqCol, generateHelm, false, chiralityEngine, false, rules);
  return ptConvertRes[0];
}


//top-menu: Bio | PolyTool | Enumerate HELM...
//name: polyToolEnumerateHelm
//description: Perform cyclization of polymers
export async function polyToolEnumerateHelmTopMenu(): Promise<void> {
  await polyToolEnumerateHelmUI(grok.shell.tv?.dataFrame.currentCell);
}

//top-menu: Bio | PolyTool | Enumerate Chem...
//name: polyToolEnumerateChem
//description: Perform cyclization of polymers
export async function polyToolEnumerateChemTopMenu(): Promise<void> {
  polyToolEnumerateChemUI();
}

//name: polyToolColumnChoice
//input: dataframe df [Input data table]
//input: column macroMolecule
export async function polyToolColumnChoice(df: DG.DataFrame, macroMolecule: DG.Column): Promise<void> {
  _setPeptideColumn(macroMolecule);
  await grok.data.detectSemanticTypes(df);
}

//name: createMonomerLibraryForPolyTool
//input: file file
export async function createMonomerLibraryForPolyTool(file: DG.FileInfo) {
  const fileContent = await file.readAsString();
  const libHandler = new PolyToolCsvLibHandler(file.fileName, fileContent);
  const libObject = await libHandler.getJson();
  const jsonFileName = file.fileName.replace(/\.csv$/, '.json');
  const jsonFileContent = JSON.stringify(libObject, null, 2);
  DG.Utils.download(jsonFileName, jsonFileContent);
}

// //name: PolyTool Converter
// //meta.icon: img/icons/structure.png
// //meta.browsePath: PolyTool
// //tags: app
// export async function ptConverterApp(): Promise<void> {
//   const view = grok.shell.v as DG.TableView;
//   const table = view.dataFrame;
//   const colNames = table.columns.names();
//   let covertableName = '';

//   for (let i = 0; i < colNames.length; i++) {
//     const col = table.columns.byName(colNames[i]);
//     if (col.semType === DG.SEMTYPE.MACROMOLECULE && col.meta.units === NOTATION.SEPARATOR) {
//       covertableName = colNames[i];
//       break;
//     }
//   }

//   if (covertableName === '')
//     grok.shell.error('To run the app open a view with convertable separator notation for macromolecules');
//   else {
//     const dialog = await getPolyToolConversionDialog();
//     dialog.show();
//   }
// }

//name: HELM Enumerator
//meta.icon: img/icons/structure.png
//meta.browsePath: Peptides | PolyTool
//tags: app
export async function ptEnumeratorHelmApp(): Promise<void> {
  await polyToolEnumerateHelmUI();
}

//name: Chem Enumerator
//meta.icon: img/icons/structure.png
//meta.browsePath: Peptides | PolyTool
//tags: app
export async function ptEnumeratorChemApp(): Promise<void> {
  polyToolEnumerateChemUI();
}


//name: Polytool Helm Enumerator dialog
//input: object cell {nullable: true}
export async function getPtHelmEnumeratorDialog(cell?: DG.Cell) {
  return polyToolEnumerateHelmUI(cell);
}

//name: Polytool Chem Enumerator dialog
//input: object cell {nullable: true}
export async function getPtChemEnumeratorDialog(cell?: DG.Cell) {
  return polyToolEnumerateChemUI(cell);
}

//name: Combine Sequences
//top-menu: Bio | PolyTool | Combine Sequences...
export async function getPolyToolCombineDialog() {
  getPTCombineDialog();
}


//name: applyNotationProviderForHarmonizedSequence
//input: column col
//input: string separator
export function applyNotationProviderForCyclized(col: DG.Column<string>, separator: string) {
  col.setTag('aligned', 'SEQ');
  col.setTag('alphabet', 'UN');
  col.setTag('.alphabetIsMultichar', 'true');
  col.meta.units = NOTATION.CUSTOM;
  col.tags[PolyToolTags.dataRole] = 'template';
  col.temp[SeqTemps.notationProvider] = new CyclizedNotationProvider(separator, _package.helmHelper);
}
