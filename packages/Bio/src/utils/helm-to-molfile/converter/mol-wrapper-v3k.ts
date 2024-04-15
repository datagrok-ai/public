import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';
import {MolfileWrapper} from './mol-wrapper';
import {RGroupHandler} from './r-group-handler';
import {MolfileAtomsV3K} from './mol-atoms-v3k';
import {MolfileBondsV3K} from './mol-bonds-v3k';

export class MolfileV3KWrapper extends MolfileWrapper {
  constructor(
    molfileV3K: string, monomerSymbol: string
  ) {
    super(monomerSymbol);
    const molfileHandler = MolfileHandler.getInstance(molfileV3K);

    this.atoms = new MolfileAtomsV3K(molfileHandler);
    this.bonds = new MolfileBondsV3K(molfileHandler);
    this.rGroups = new RGroupHandler(molfileHandler, this.atoms, this.bonds);

    this.shiftMonomerToDefaultPosition();
  }
}

