import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {OligoTranslatorUI, OligoPatternUI, OligoStructureUI, AppUI, AppMultiView} from './view/view';
import {tryCatch} from './model/helpers';
import {LIB_PATH, DEFAULT_LIB_FILENAME} from './model/data-loading-utils/const';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {getJsonData} from './model/data-loading-utils/json-loader';
import {SequenceToMolfileConverter} from './model/sequence-to-structure-utils/sequence-to-molfile';
import {linkStrandsV3000} from './model/sequence-to-structure-utils/mol-transformations';
import {MonomerLibWrapper} from './model/monomer-lib/lib-wrapper';
import {FormatDetector} from './model/parsing-validation/format-detector';
import {SequenceValidator} from './model/parsing-validation/sequence-validator';
import {demoOligoTranslatorUI, demoOligoPatternUI} from './demo/demo-st-ui';
import {FormatConverter} from './model/format-translation/format-converter';
import {COMBINED_APP_NAME, PATTERN_APP_NAME, STRUCTRE_APP_NAME, TRANSLATOR_APP_NAME} from './view/const/view';

class StPackage extends DG.Package {
  private _monomerLib?: IMonomerLib;

  get monomerLib(): IMonomerLib {
    if (!this._monomerLib)
      throw new Error ('Monomer lib not loaded')
    return this._monomerLib!;
  }

  public async initMonomerLib(): Promise<void> {
    if (this._monomerLib !== undefined)
      return;

    const pi: DG.TaskBarProgressIndicator = DG.TaskBarProgressIndicator.create(
      `Initializing ${COMBINED_APP_NAME} monomer library ...`);
    await tryCatch(async () => {
      const libHelper: IMonomerLibHelper = await getMonomerLibHelper();
      this._monomerLib = await libHelper.readLibrary(LIB_PATH, DEFAULT_LIB_FILENAME);
    }, () => pi.close());
  }
}

export const _package: StPackage = new StPackage();

//name: Oligo Toolkit
//tags: app
export async function oligoToolkitApp(): Promise<void> {
  const pi: DG.TaskBarProgressIndicator = DG.TaskBarProgressIndicator.create(`Loading ${COMBINED_APP_NAME}...`);

  await tryCatch(async () => {
    await initSequenceTranslatorLibData();
    const multiView = new AppMultiView();
    multiView.createLayout();
  }, () => pi.close());
}

async function createAppLayout(appUI: AppUI, appName: string): Promise<void> {
  const pi: DG.TaskBarProgressIndicator = DG.TaskBarProgressIndicator.create(`Loading ${appName}...`);

  await tryCatch(async () => {
    await initSequenceTranslatorLibData();
    await appUI.createLayout();
  }, () => pi.close());
}

//name: Oligo Translator
//tags: app
export async function oligoTranslatorApp(): Promise<void> {
  await initSequenceTranslatorLibData();
  const appUI = new OligoTranslatorUI(DG.View.create());
  await createAppLayout(appUI, TRANSLATOR_APP_NAME);
}

//name: Oligo Pattern
//tags: app
export async function oligoPatternApp(): Promise<void> {
  await initSequenceTranslatorLibData();
  const appUI = new OligoPatternUI(DG.View.create());
  createAppLayout(appUI, PATTERN_APP_NAME);
}

//name: Oligo Structure
//tags: app
export async function oligoStructureApp(): Promise<void> {
  await initSequenceTranslatorLibData();
  const appUI = new OligoStructureUI(DG.View.create());
  createAppLayout(appUI, STRUCTRE_APP_NAME);
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

//name: demoOligoTranslator
//meta.demoPath: Bioinformatics | Oligo Toolkit | Translator
//description: Translate oligonucleotide sequences across various formats accepted by different synthesizers
//meta.path: /apps/Tutorials/Demo/Bioinformatics/Oligonucleotide%20Sequence:%20Translate
export async function demoTranslateSequence(): Promise<void> {
  await demoOligoTranslatorUI();
}

//name: demoOligoPattern
//meta.demoPath: Bioinformatics | Oligo Toolkit | Pattern
//description: Design a modification pattern for an oligonucleotide sequence
//meta.path:%20/apps/Tutorials/Demo/Bioinformatics/Oligonucleotide%20Sequence:%20Visualize%20duplex
export async function demoDesignPattern(): Promise<void> {
  await demoOligoPatternUI();
}

//name: demoVisualizeDuplex
//meta.demoPath: Bioinformatics | Oligonucleotide Sequence: Visualize duplex
//description: Visualize duplex and save SDF
//meta.path:%20/apps/Tutorials/Demo/Bioinformatics/Oligonucleotide%20Sequence:%20Visualize%20duplex
export async function demoVisualizeDuplex(): Promise<void> {
  // await demoVisualizeDuplexUI();
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
