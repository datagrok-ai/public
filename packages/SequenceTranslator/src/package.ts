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

export * from './package.g';
export const _package: OligoToolkitPackage = new OligoToolkitPackage({debug: true}/**/);

let initSequenceTranslatorPromise: Promise<void> | null = null;

async function initSequenceTranslatorInt(): Promise<void> {
  const [helmHelper] = await Promise.all([
    getHelmHelper(),
  ]);
  _package.completeInit(helmHelper);
}

export class PackageFunctions{
  @grok.decorators.app({
    icon: 'img/icons/toolkit.png',
    browsePath: 'Peptides | Oligo Toolkit',
    name: 'Oligo Toolkit'
  })
  static async oligoToolkitApp(): Promise<DG.ViewBase> {
  
    await _package.initLibData();
    const externalViewFactories = await getExternalAppViewFactories(_package);
    if (!externalViewFactories)
      throw new Error('External app view factories not loaded');
    const appUI = new CombinedAppUI(externalViewFactories!, _package);
    const view = await appUI.getAppView();
    return view;
  }

    
  @grok.decorators.init()
  static async init(): Promise<void> {
  
    if (initSequenceTranslatorPromise === null)
      _package.startInit(initSequenceTranslatorPromise = initSequenceTranslatorInt());

    return initSequenceTranslatorPromise;
  }

  @grok.decorators.app({
    icon: 'img/icons/translator.png',
    browsePath: 'Peptides | Oligo Toolkit',
    name: 'Oligo Translator'
  })
  static async oligoTranslatorApp(): Promise<DG.ViewBase> {
  
    const view = await getSpecifiedAppView(APP_NAME.TRANSLATOR);
    return view;
  }


  @grok.decorators.app({
    icon: 'img/icons/pattern.png',
    browsePath: 'Peptides | Oligo Toolkit',
    name: 'Oligo Pattern'
  })
  static async oligoPatternApp(): Promise<DG.ViewBase> {
  
    const view = await getSpecifiedAppView(APP_NAME.PATTERN);
    return view;
  }


  @grok.decorators.app({
    icon: 'img/icons/structure.png',
    browsePath: 'Peptides | Oligo Toolkit',
    name: 'Oligo Structure'
  })
  static async oligoStructureApp(): Promise<DG.ViewBase> {
  
    const view = await getSpecifiedAppView(APP_NAME.STRUCTURE);
    return view;
  }


  @grok.decorators.func({outputs: [{type: 'object', name: 'result'}]})
  static async getTranslationHelper(): Promise<ITranslationHelper> {
  
    await _package.initLibData();
    return _package;
  }


  @grok.decorators.func({outputs: [{type: 'object', name: 'result'}]})
  static getCodeToWeightsMap(): { [key: string]: number } {
  
    const monomerLibWrapper = _package.monomerLibWrapper;
    const map = monomerLibWrapper.getCodesToWeightsMap();
    return Object.fromEntries(map);
  }


  @grok.decorators.func()
  static validateSequence(
    sequence: string): boolean {
  
    const validator = _package.createSequenceValidator(sequence);
    const format = _package.createFormatDetector(sequence).getFormat();
    return (format === null) ? false : validator.isValidSequence(format!);
  }


  @grok.decorators.func({
    name: 'validateSequence'
  })
  static getMolfileFromGcrsSequence(
    sequence: string,
    invert: boolean): string {
  
    return (new SequenceToMolfileConverter(sequence, invert, 'GCRS')).convert();
  }


  @grok.decorators.func()
  static linkStrands(
    @grok.decorators.param({type:'object'})  strands: { senseStrands: string[], antiStrands: string[] }): string {
    return linkStrandsV3000(strands, true);
  }


  @grok.decorators.func({
    meta: {
      demoPath: 'Bioinformatics | Oligo Toolkit | Translator',
      path: '/apps/Tutorials/Demo/Bioinformatics/Oligonucleotide%20Sequence:%20Translate',
      demoSkip: 'GROK-14320'
    },
    name: 'demoOligoTranslator',
    description: 'Translate oligonucleotide sequences across various formats accepted by different synthesizers'
  })
  static async demoTranslateSequence(): Promise<void> {
    await demoOligoTranslatorUI();
  }


  @grok.decorators.func({
    meta: {
      demoPath: 'Bioinformatics | Oligo Toolkit | Pattern',
      path: '%20/apps/Tutorials/Demo/Bioinformatics/Oligonucleotide%20Sequence:%20Visualize%20duplex'
    },
    description: 'Design a modification pattern for an oligonucleotide sequence'
  })
  static async demoOligoPattern(): Promise<void> {
  
    await demoOligoPatternUI();
  }


  @grok.decorators.func({
    meta: {
      demoPath: 'Bioinformatics | Oligo Toolkit | Structure',
      path: '%20/apps/Tutorials/Demo/Bioinformatics/Oligonucleotide%20Sequence:%20Visualize%20duplex'
    },
    description: 'Visualize duplex and save SDF'
  })
  static async demoOligoStructure(): Promise<void> {
  
    await demoOligoStructureUI();
  }


  @grok.decorators.func()
  static async translateOligonucleotideSequence(
    sequence: string, sourceFormat: string, targetFormat: string
  ): Promise<string> {
  
    await _package.initLibData();
    return _package.createFormatConverter(sequence, sourceFormat).convertTo(targetFormat);
  }

  @grok.decorators.func({
    'top-menu': 'Bio | PolyTool | Convert...',
    name: 'polyToolConvert',
    description: 'Perform cyclization of polymers'
  })
  static async polyToolConvertTopMenu(): Promise<void> {
  
    await polyToolConvertUI();
  }

  @grok.decorators.editor()
  static async getPolyToolConvertEditor(
    call: DG.FuncCall): Promise<DG.Column<string> | null> {
  
    const funcEditor = await PolyToolConvertFuncEditor.create(call);
    return await funcEditor.showDialog();
  }


  @grok.decorators.func({
    editor: 'SequenceTranslator:getPolyToolConvertEditor', 
  })
  static async polyToolConvert2(
    table: DG.DataFrame,
    @grok.decorators.param({options: { caption: 'Sequence'}}) seqCol: DG.Column, 
    @grok.decorators.param({options: { initialValue: 'true'}}) generateHelm: boolean, 
    @grok.decorators.param({options: { initialValue: 'true'}}) chiralityEngine: boolean, 
    @grok.decorators.param({type: 'object'}) rules: string[]
  ): Promise<DG.Column<string>> {
  
    const ptConvertRes = await polyToolConvert(seqCol, generateHelm, false, chiralityEngine, false, rules);
    return ptConvertRes[0];
  }



  @grok.decorators.func({
    'top-menu': 'Bio | PolyTool | Enumerate HELM...',
    name: 'polyToolEnumerateHelm',
    description: 'Perform cyclization of polymers'
  })
  static async polyToolEnumerateHelmTopMenu(): Promise<void> {
  
    await polyToolEnumerateHelmUI(grok.shell.tv?.dataFrame.currentCell);
  }


  @grok.decorators.func({
    'top-menu': 'Bio | PolyTool | Enumerate Chem...',
    name: 'polyToolEnumerateChem',
    description: 'Perform cyclization of polymers'
  })
  static async polyToolEnumerateChemTopMenu(): Promise<void> {
  
    polyToolEnumerateChemUI();
  }


  @grok.decorators.func()
  static async polyToolColumnChoice(
    @grok.decorators.param({options: {description: 'Input data table'}}) df: DG.DataFrame,
    macroMolecule: DG.Column): Promise<void> {
  
    _setPeptideColumn(macroMolecule);
    await grok.data.detectSemanticTypes(df);
  }


  @grok.decorators.func()
  static async createMonomerLibraryForPolyTool(
    file: DG.FileInfo) {
  
    const fileContent = await file.readAsString();
    const libHandler = new PolyToolCsvLibHandler(file.fileName, fileContent);
    const libObject = await libHandler.getJson();
    const jsonFileName = file.fileName.replace(/\.csv$/, '.json');
    const jsonFileContent = JSON.stringify(libObject, null, 2);
    DG.Utils.download(jsonFileName, jsonFileContent);
  }


  @grok.decorators.func({
    meta: {
      icon: 'img/icons/structure.png',
      browsePath: 'Peptides | PolyTool'
    },
    tags: ['app'],
    name: 'HELM Enumerator'
  })
  static async ptEnumeratorHelmApp(): Promise<void> {
  
    await polyToolEnumerateHelmUI();
  }


  @grok.decorators.func({
    meta: {
      icon: 'img/icons/structure.png',
      browsePath: 'Peptides | PolyTool'
    },
    tags: ['app'],
    name: 'Chem Enumerator'
  })
  static async ptEnumeratorChemApp(): Promise<void> {
  
    polyToolEnumerateChemUI();
  }



  @grok.decorators.func({
    name: 'Polytool Helm Enumerator dialog'
  })
  static async getPtHelmEnumeratorDialog(
    @grok.decorators.param({type:'object',options:{nullable:true}})  cell?: DG.Cell) {
  
    return polyToolEnumerateHelmUI(cell);
  }


  @grok.decorators.func({
    name: 'Polytool Chem Enumerator dialog'
  })
  static async getPtChemEnumeratorDialog(
    @grok.decorators.param({type:'object',options:{nullable:true}})  cell?: DG.Cell) {
  
    return polyToolEnumerateChemUI(cell);
  }


  @grok.decorators.func({
    name: 'Combine Sequences',
    'top-menu': 'Bio | PolyTool | Combine Sequences...'
  })
  static async getPolyToolCombineDialog() {
  
    getPTCombineDialog();
  }



  @grok.decorators.func({
    name: 'applyNotationProviderForHarmonizedSequence'
  })
  static applyNotationProviderForCyclized(
    @grok.decorators.param({type:'column'}) col: DG.Column<string>,
    separator: string) {
  
    col.setTag('aligned', 'SEQ');
    col.setTag('alphabet', 'UN');
    col.setTag('.alphabetIsMultichar', 'true');
    col.meta.units = NOTATION.CUSTOM;
    col.tags[PolyToolTags.dataRole] = 'template';
    col.temp[SeqTemps.notationProvider] = new CyclizedNotationProvider(separator, _package.helmHelper);
  }
}



async function getSpecifiedAppView(appName: string): Promise<DG.ViewBase> {
  await _package.initLibData();
  const appUI = getSpecifiedAppUI(appName, _package);
  const view = await appUI.getAppView();
  return view;
}