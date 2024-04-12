import {_package} from '../../package';
import {R_GROUP_ELEMENT_SYMBOL} from './const';

export class MolfileAtoms {
  constructor(atomLines: string[]) {
    this.rawAtomLines = atomLines;
    this.coordinates = this.rawAtomLines.map((line: string) => {
      const x = parseFloat(line.substring(0, 10));
      const y = parseFloat(line.substring(10, 20));
      return {x, y};
    });
  }

  private coordinates: {x: number, y: number}[] = [];
  private rawAtomLines: string[] = [];

  get atomCoordinates(): {x: number, y: number}[] {
    return this.coordinates;
  }

  get atomLines(): string[] {
    return this.rawAtomLines.map((line: string, idx: number) => {
      const coordinates = this.coordinates[idx];
      const x = coordinates.x.toFixed(4).padStart(10, ' ');
      const y = coordinates.y.toFixed(4).padStart(10, ' ');
      return `${x}${y}${line.substring(20)}`;
    });
  }

  replaceElementSymbol(atomIdx: number, newElementSymbol: string): void {
    this.rawAtomLines[atomIdx] = this.rawAtomLines[atomIdx].replace(R_GROUP_ELEMENT_SYMBOL, newElementSymbol);
  }

  deleteAtoms(indices: number[]): void {
    this.coordinates = this.coordinates.filter((_, idx) => !indices.includes(idx));
    this.rawAtomLines = this.rawAtomLines.filter((_, idx) => !indices.includes(idx));
  }

  shift(shift: {x: number, y: number}): void {
    this.coordinates = this.coordinates.map((coordinates) => {
      const newX = coordinates.x + shift.x;
      const newY = coordinates.y + shift.y;
      if (isNaN(newX) || isNaN(newY))
        throw new Error(`Cannot shift coordinates by ${shift.x}, ${shift.y}`);
      return {x: newX, y: newY};
    });
  }

  rotate(angle: number): void {
    this.coordinates = this.coordinates.map((coordinates) => {
      const x = coordinates.x;
      const y = coordinates.y;
      const newX = x * Math.cos(angle) - y * Math.sin(angle);
      const newY = x * Math.sin(angle) + y * Math.cos(angle);
      if (isNaN(newX) || isNaN(newY))
        throw new Error(`Cannot rotate coordinates by ${angle}`);
      return {x: newX, y: newY};
    });
  }

  getRGroupAtomicIndices(): number[] {
    return this.rawAtomLines.map((line: string, idx: number) => {
      if (line.includes(R_GROUP_ELEMENT_SYMBOL))
        return idx;
    }).filter((idx) => idx !== undefined) as number[];
  }
}

