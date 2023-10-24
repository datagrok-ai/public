import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {OligoTranslatorUI, OligoPatternUI, OligoStructureUI, CompundAppUI} from './view/app-ui';
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
import {demoOligoTranslatorUI, demoOligoPatternUI, demoOligoStructureUI} from './demo/demo-st-ui';
import {FormatConverter} from './model/format-translation/format-converter';
import {APP} from './view/const/view';
import {getExternalAppViewFactories} from './plugins/mermade';

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
      `Initializing ${APP.COMBINED} monomer library ...`);
    await tryCatch(async () => {
      const libHelper: IMonomerLibHelper = await getMonomerLibHelper();
      this._monomerLib = await libHelper.readLibrary(LIB_PATH, DEFAULT_LIB_FILENAME);
    }, () => pi.close());
  }
}

export const _package: StPackage = new StPackage();

//name: Oligo Toolkit
//meta.icon: img/icons/toolkit.png
//tags: app
export async function oligoToolkitApp(): Promise<void> {
  await initSequenceTranslatorLibData();
  const externalViewFactories = await getExternalAppViewFactories();
  if (!externalViewFactories)
    throw new Error('External app view factories not loaded');
  const appUI = new CompundAppUI(externalViewFactories!);
  await appUI.createAppLayout();
}

//name: Oligo Translator
//meta.icon: img/icons/translator.png
//tags: app
export async function oligoTranslatorApp(): Promise<void> {
  await initSequenceTranslatorLibData();
  const appUI = new OligoTranslatorUI();
  await appUI.createAppLayout();
}

//name: Oligo Pattern
//meta.icon: img/icons/pattern.png
//tags: app
export async function oligoPatternApp(): Promise<void> {
  await initSequenceTranslatorLibData();
  const appUI = new OligoPatternUI();
  await appUI.createAppLayout();
}

//name: Oligo Structure
//meta.icon: img/icons/structure.png
//tags: app
export async function oligoStructureApp(): Promise<void> {
  await initSequenceTranslatorLibData();
  const appUI = new OligoStructureUI();
  await appUI.createAppLayout();
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
export async function translateOligonucleotideSequence(sequence: string, sourceFormat: string, targetFormat: string): Promise<string> {
  await initSequenceTranslatorLibData();
  return (new FormatConverter(sequence, sourceFormat)).convertTo(targetFormat);
}
