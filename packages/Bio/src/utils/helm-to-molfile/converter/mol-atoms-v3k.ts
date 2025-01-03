import {MolfileHandlerBase} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler-base';
import {MolfileAtoms} from './mol-atoms';

const PRECISION = 4;

export class MolfileAtomsV3K extends MolfileAtoms {
  constructor(private molfileHandler: MolfileHandlerBase) {
    super();
    this.rawAtomLines = molfileHandler.getAtomLines();
    this.coordinates = this.getCoordinates();
  }

  private getCoordinates(): {x: number, y: number}[] {
    const x = this.molfileHandler.x;
    const y = this.molfileHandler.y;
    return Array.from(x).map((xCoord, idx) => {
      return {x: xCoord, y: y[idx]};
    });
  }


  get atomLines(): string[] {
    // todo: optimize, optionally port to molfile-handler
    const coordinateRegex = /^(M  V30 [^-]*)(-?\d+\.\d+)( )(-?\d+\.\d+)( -?\d+\.\d+.*)$/;
    const rGroupsRegex = /\sRGROUPS=\(\d+(\s+\d+)*\)/;

    return this.rawAtomLines.map((line: string, idx: number) => {
      const coordinates = this.coordinates[idx];
      const x = coordinates.x.toFixed(PRECISION) + '00';
      const y = coordinates.y.toFixed(PRECISION) + '00';

      return line.replace(coordinateRegex, (match, p1, p2, p3, p4, p5) => {
        return p1 + x + p3 + y + p5;
      }).replace(rGroupsRegex, '');
    });
  }

  replaceRGroupSymbolByElement(atomIdx: number, newElementSymbol: string): void {
    super.replaceRGroupSymbolByElement(atomIdx, newElementSymbol);
    // rdkit can generate (out of thin air) masses for r groups, so we need to remove them as well.
    //they are at the end of the line after coordinates and other data
    const lineInfo = this.rawAtomLines[atomIdx].substring(3).split(' ');
    if (lineInfo.length > 7)
      this.rawAtomLines[atomIdx] = `M  ${lineInfo.slice(0, 7).join(' ')}`;
  }
}

