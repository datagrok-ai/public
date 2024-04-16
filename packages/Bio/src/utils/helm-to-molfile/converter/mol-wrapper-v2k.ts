import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';
import {MolfileWrapper} from './mol-wrapper';
import {RGroupHandler} from './r-group-handler';
import {MolfileAtomsV2K} from './mol-atoms-v2k';
import {MolfileBondsV2K} from './mol-bonds-v2k';

export class MolfileV2KWrapper extends MolfileWrapper {
  constructor(
    molfileV2K: string, monomerSymbol: string
  ) {
    super(monomerSymbol);
    const molfileHandler = MolfileHandler.getInstance(molfileV2K);

    this.atoms = new MolfileAtomsV2K(molfileHandler);
    this.bonds = new MolfileBondsV2K(molfileHandler);
    this.rGroups = new RGroupHandler(molfileHandler, this.atoms, this.bonds);

    this.shiftMonomerToDefaultPosition();
  }
}

