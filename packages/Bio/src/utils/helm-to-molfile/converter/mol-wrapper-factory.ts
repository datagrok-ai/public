import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';
import {MolfileWrapper} from './mol-wrapper';
import {MolfileV2KWrapper} from './mol-wrapper-v2k';
import {MolfileV3KWrapper} from './mol-wrapper-v3k';

export class MolfileWrapperFactory {
  static getInstance(molfile: string, monomerSymbol: string): MolfileWrapper {
    if (MolfileHandler.isMolfileV2K(molfile))
      return new MolfileV2KWrapper(molfile, monomerSymbol) as MolfileWrapper;
    else if (MolfileHandler.isMolfileV3K(molfile))
      return new MolfileV3KWrapper(molfile, monomerSymbol) as MolfileWrapper;
    else
      throw new Error('Unsupported molfile version');
  }
}

