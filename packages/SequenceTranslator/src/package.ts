import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {SequenceTranslatorUI} from './view/view';
import {LIB_PATH, DEFAULT_LIB_FILENAME} from './model/data-loading-utils/const';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {getJsonData} from './model/data-loading-utils/json-loader';
import {SequenceToMolfileConverter} from './model/sequence-to-structure-utils/sequence-to-molfile';
import {linkStrandsV3000} from './model/sequence-to-structure-utils/mol-transformations';
import {MonomerLibWrapper} from './model/monomer-lib/lib-wrapper';
import {FormatDetector} from './model/parsing-validation/format-detector';
import {SequenceValidator} from './model/parsing-validation/sequence-validator';
import {demoDesignPatternUI, demoVisualizeDuplexUI, demoTranslateSequenceUI} from './demo/demo-st-ui';
import {FormatConverter} from './model/format-translation/format-converter';

class StPackage extends DG.Package {
  private _monomerLib?: IMonomerLib;

  get monomerLib(): IMonomerLib {
    if (!this._monomerLib)
      throw new Error ('ST: monomer lib not loaded')
    return this._monomerLib!;
  }

  public async initMonomerLib(): Promise<void> {
    if (this._monomerLib !== undefined)
      return;

    const pi: DG.TaskBarProgressIndicator = DG.TaskBarProgressIndicator.create(
      'Initializing Sequence Translator monomer library ...');
    try {
      const libHelper: IMonomerLibHelper = await getMonomerLibHelper();
      this._monomerLib = await libHelper.readLibrary(LIB_PATH, DEFAULT_LIB_FILENAME);
    } catch (err: any) {
      const errMsg: string = err.hasOwnProperty('message') ? err.message : err.toString();
      throw new Error('Sequence Translator: Loading monomer library error: ' + errMsg);
    } finally {
      pi.close();
    }
  }
}

export const _package: StPackage = new StPackage();

//name: Sequence Translator
//tags: app
export async function sequenceTranslatorApp(): Promise<void> {
  const pi: DG.TaskBarProgressIndicator = DG.TaskBarProgressIndicator.create('Loading Sequence Translator app ...');
  let currentView = grok.shell.v.root;
  ui.setUpdateIndicator(currentView,true);
  try {
    await initSequenceTranslatorLibData();
    const v = new SequenceTranslatorUI();
    await v.createLayout();
  } catch (err: any) {
    const errMsg: string = err.hasOwnProperty('message') ? err.message : err.toString();
    grok.shell.error(`Loading Sequence Translator application error: ` + errMsg);
    throw err;
  } finally {
    pi.close();
  }
  ui.setUpdateIndicator(currentView,false);
}

//name: initSequenceTranslatorLibData
export async function initSequenceTranslatorLibData(): Promise<void> {
  await getJsonData();
  await _package.initMonomerLib();
}

//name: getCodeToWeightsMap
//output: object result
export function getCodeToWeightsMap(): {[key: string]: number} {
  const map = MonomerLibWrapper.getInstance().getCodesToWeightsMap();
  return Object.fromEntries(map);
}

//name: validateSequence
//input: string sequence
//output: bool result
export function validateSequence(sequence: string): boolean {
  const validator = new SequenceValidator(sequence);
  const format = (new FormatDetector(sequence).getFormat());
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

//name: demoTranslateSequence
//meta.demoPath: Bioinformatics | Oligonucleotide Sequence: Translate
//description: Translate oligonucleotide sequences across various formats accepted by different synthesizers
//meta.path: /apps/Tutorials/Demo/Bioinformatics/Oligonucleotide%20Sequence:%20Translate
export async function demoTranslateSequence(): Promise<void> {
  await demoTranslateSequenceUI();
}

//name: demoDesignPattern
//meta.demoPath: Bioinformatics | Oligonucleotide Sequence: Design
//description: Design a modification pattern for an oligonucleotide sequence
//meta.path:%20/apps/Tutorials/Demo/Bioinformatics/Oligonucleotide%20Sequence:%20Visualize%20duplex
export async function demoDesignPattern(): Promise<void> {
  await demoDesignPatternUI();
}

//name: demoVisualizeDuplex
//meta.demoPath: Bioinformatics | Oligonucleotide Sequence: Visualize duplex
//description: Visualize duplex and save SDF
//meta.path:%20/apps/Tutorials/Demo/Bioinformatics/Oligonucleotide%20Sequence:%20Visualize%20duplex
export async function demoVisualizeDuplex(): Promise<void> {
  await demoVisualizeDuplexUI();
}

//name: translateOligonucleotideSequence
//input: string sequence
//input: string sourceFormat
//input: string targetFormat
//output: string result
export async function translateOligonucleotideSequence(sequence: string, sourceFormat: string, targetFormat: string): Promise<string> {
  await initSequenceTranslatorLibData();
  return (new FormatConverter(sequence, sourceFormat)).convertTo(targetFormat);
}
