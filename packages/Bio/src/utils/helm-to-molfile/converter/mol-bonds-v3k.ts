import {MolfileHandlerBase} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler-base';
import {MolfileBonds} from './mol-bonds';

export class MolfileBondsV3K extends MolfileBonds {
  constructor(private molfileHandler: MolfileHandlerBase) {
    super();
    this.rawBondLines = molfileHandler.getBondLines();
    this.bondedAtomPairs = this.getBondedAtomPairs();
  }

  private getBondedAtomPairs(): number[][] {
    const bondedAtoms = this.molfileHandler.pairsOfBondedAtoms;
    return bondedAtoms.map((pair) => [pair[0], pair[1]]);
  }

  /** Get bond lines with new values for bonded atoms  */
  getBondLines(): string[] {
    // todo: optimize
    const regex = /^(M\s+V30\s+\d+\s+\d+\s+)(\d+)(\s+)(\d+)(.*)$/;
    return this.bondedAtomPairs.map((bondedPair, idx) => {
      if (bondedPair.some((atom) => atom === -1))
        throw new Error(`Bonded pair ${bondedPair} contains -1`);
      const result = this.rawBondLines[idx].replace(regex, (match, p1, p2, p3, p4, p5) => {
        return `${p1}${bondedPair[0]}${p3}${bondedPair[1]}${p5}`;
      });
      return result;
    });
  }
}

