import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';
import {MolfileWrapperBase} from './mol-wrapper-base';
import {MolfileV2KWrapper} from './mol-v2k-wrapper';

export class MolfileWrapper {
  static getInstance(molfile: string, monomerSymbol: string): MolfileWrapperBase {
    if (MolfileHandler.isMolfileV2K(molfile))
      return new MolfileV2KWrapper(molfile, monomerSymbol) as MolfileWrapperBase;
    else
      throw new Error('Unsupported molfile version');
  }
}

