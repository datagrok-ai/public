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
import {demoSequenceTranslatorUI} from './demo/demo-st-ui';

class StPackage extends DG.Package {
  private _monomerLib?: IMonomerLib;

  get monomerLib(): IMonomerLib {
    if (!this._monomerLib)
      throw new Error ('ST: monomer lib not loaded')
    return this._monomerLib!;
  }

  public async initMonomerLib(): Promise<void> {
    if (this._monomerLib === undefined) {
      const pi: DG.TaskBarProgressIndicator = DG.TaskBarProgressIndicator.create(
        'Initializing Sequence Translator monomer library ...');
      try {
        const libHelper: IMonomerLibHelper = await getMonomerLibHelper();
        this._monomerLib = await libHelper.readLibrary(LIB_PATH, DEFAULT_LIB_FILENAME);
      } catch (err: any) {
        const errMsg: string = err.hasOwnProperty('message') ? err.message : err.toString();
        throw new Error('ST: Loading monomer library error: ' + errMsg);
      } finally {
        pi.close();
      }
    }
  }
}

export const _package: StPackage = new StPackage();

//name: Sequence Translator
//tags: app
export async function sequenceTranslatorApp(): Promise<void> {
  const pi: DG.TaskBarProgressIndicator = DG.TaskBarProgressIndicator.create('Loading Sequence Translator app ...');

  try {
    await getJsonData();
    await _package.initMonomerLib();
    const v = new SequenceTranslatorUI();
    await v.createLayout();
  } catch (err: any) {
    const errMsg: string = err.hasOwnProperty('message') ? err.message : err.toString();
    grok.shell.error(`Loading Sequence Translator application error: ` + errMsg);
    throw err;
  } finally {
    pi.close();
  }
}

//name: getCodeToNameMap
export function getCodeToNameMap(): Map<string, number> {
  return MonomerLibWrapper.getInstance().getCodesToWeightsMap();
}

//name: validateSequence
//input: string sequence
export function validateSequence(sequence: string): boolean {
  const validator = new SequenceValidator(sequence);
  const format = (new FormatDetector(sequence).getFormat());
  return (format === null) ? false : validator.isValidSequence(format!);
}

//name: validateSequence
//input: string sequence
//input: bool invert
export function getMolfileFromGcrsSequence(sequence: string, invert: boolean): string {
  return (new SequenceToMolfileConverter(sequence, invert, 'GCRS')).convert();
}

//name: linkStrands
//input: object strands
//input: bool invert
export function linkStrands(strands: { senseStrands: string[], antiStrands: string[] }): string {
  return linkStrandsV3000(strands, true);
}

// demoSequenceTranslator
//name: demoSequenceTranslator
//meta.demoPath: Bioinformatics | Sequence Design
//description: Sequence Translator is an application for design and visualization of oligonucleotide sequences
//meta.path: /apps/Tutorials/Demo/Bioinformatics/Sequence,%20Translator
//meta.isDemoScript: True
export async function demoSequenceTranslator(): Promise<void> {
  await demoSequenceTranslatorUI();
}
