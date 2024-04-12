import {PositionInBonds} from './types';

export class MolfileBonds {
  constructor(bondLines: string[]) {
    this.rawBondLines = bondLines;
    this.bondedPairs = this.rawBondLines.map((line: string) => {
      const firstAtom = parseInt(line.substring(0, 3));
      const secondAtom = parseInt(line.substring(3, 6));
      return [firstAtom, secondAtom];
    });
  }

  private bondedPairs: number[][] = [];
  private rawBondLines: string[] = [];

  /** Get bond lines with new values for bonded atoms  */
  getBondLines(): string[] {
    return this.bondedPairs.map((bondedPair, idx) => {
      if (bondedPair.some((atom) => atom === -1))
        throw new Error(`Bonded pair ${bondedPair} contains -1`);
      return `${bondedPair[0].toString().padStart(3, ' ')}${
        bondedPair[1].toString().padStart(3, ' ')
      }${this.rawBondLines[idx].substring(6)}`;
    });
  }

  get bondedAtoms(): number[][] {
    return this.bondedPairs;
  }

  deleteBondLines(indices: number[]): void {
    this.rawBondLines = this.rawBondLines.filter((_, idx) => !indices.includes(idx));
    this.bondedPairs = this.bondedPairs.filter((_, idx) => !indices.includes(idx));
  }

  /** Atom id starts from 1  */
  getPositionsInBonds(atomId: number): PositionInBonds[] {
    const positions: PositionInBonds[] = [];
    this.bondedPairs.forEach((bondedPair, bondLineIdx) => {
      bondedPair.forEach((atom, nodeIdx) => {
        if (atom === atomId)
          positions.push({bondLineIdx, nodeIdx});
      });
    });
    return positions;
  }

  replacePositionsInBondsByDummy(positions: PositionInBonds[], dummy?: number): void {
    if (dummy === undefined)
      dummy = -1;
    positions.forEach((position) => {
      const {bondLineIdx, nodeIdx} = position;
      this.bondedPairs[bondLineIdx][nodeIdx] = dummy!;
    });
  }

  removeAtomIdFromBonds(atomId: number): void {
    this.bondedPairs = this.bondedPairs.map((bondedPair) => {
      return bondedPair.map((id) => {
        if (id > atomId)
          return id - 1;
        return id;
      });
    });
  }

  shift(shift: number): void {
    this.bondedPairs = this.bondedPairs.map((bondedPair) => {
      return bondedPair.map((id) => id + shift);
    });
  }
}

