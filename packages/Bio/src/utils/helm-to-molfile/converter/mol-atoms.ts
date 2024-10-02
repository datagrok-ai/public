import {R_GROUP_ELEMENT_SYMBOL} from '@datagrok-libraries/chem-meta/src/formats/molfile-const';

export abstract class MolfileAtoms {
  protected coordinates: {x: number, y: number}[] = [];
  protected rawAtomLines: string[] = [];

  get count(): number { return this.coordinates.length; }

  get atomCoordinates(): {x: number, y: number}[] {
    return this.coordinates;
  }

  abstract get atomLines(): string[];

  replaceRGroupSymbolByElement(atomIdx: number, newElementSymbol: string): void {
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
}

