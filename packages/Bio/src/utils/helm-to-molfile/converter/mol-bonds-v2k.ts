import {MolfileHandlerBase} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler-base';
import {MolfileBonds} from './mol-bonds';

export class MolfileBondsV2K extends MolfileBonds {
  constructor(molfileHandler: MolfileHandlerBase) {
    super();
    this.rawBondLines = molfileHandler.getBondLines();
    this.bondedAtomPairs = this.rawBondLines.map((line: string) => {
      const firstAtom = parseInt(line.substring(0, 3));
      const secondAtom = parseInt(line.substring(3, 6));
      return [firstAtom, secondAtom];
    });
  }

  /** Get bond lines with new values for bonded atoms  */
  getBondLines(): string[] {
    return this.bondedAtomPairs.map((bondedPair, idx) => {
      if (bondedPair.some((atom) => atom === -1))
        throw new Error(`Bonded pair ${bondedPair} contains -1`);
      return `${bondedPair[0].toString().padStart(3, ' ')}${
        bondedPair[1].toString().padStart(3, ' ')
      }${this.rawBondLines[idx].substring(6)}`;
    });
  }
}

