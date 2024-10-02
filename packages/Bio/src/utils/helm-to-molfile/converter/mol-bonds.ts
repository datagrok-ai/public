import {PositionInBonds} from './types';

export abstract class MolfileBonds {
  protected bondedAtomPairs: number[][] = [];
  protected rawBondLines: string[] = [];

  public get count(): number { return this.bondedAtomPairs.length;}

  /** Get bond lines with new values for bonded atoms  */
  abstract getBondLines(): string[];

  get bondedAtoms(): number[][] {
    return this.bondedAtomPairs;
  }

  deleteBondLines(indices: number[]): void {
    this.rawBondLines = this.rawBondLines.filter((_, idx) => !indices.includes(idx));
    this.bondedAtomPairs = this.bondedAtomPairs.filter((_, idx) => !indices.includes(idx));
  }

  /** Atom id starts from 1  */
  getPositionsInBonds(atomId: number): PositionInBonds[] {
    const positions: PositionInBonds[] = [];
    this.bondedAtomPairs.forEach((bondedPair, bondLineIdx) => {
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
      this.bondedAtomPairs[bondLineIdx][nodeIdx] = dummy!;
    });
  }

  removeAtomIdFromBonds(atomId: number): void {
    this.bondedAtomPairs = this.bondedAtomPairs.map((bondedPair) => {
      return bondedPair.map((id) => {
        if (id > atomId)
          return id - 1;
        return id;
      });
    });
  }

  shift(shift: number): void {
    this.bondedAtomPairs = this.bondedAtomPairs.map((bondedPair) => {
      return bondedPair.map((id) => id + shift);
    });
  }
}

