import {MolfileHandlerBase} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler-base';
import {MolfileAtoms} from './mol-atoms';

const PRECISION = 6;

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
    // todo: optimize, port to molfile-handler
    const coordinateRegex = /^(M  V30 .*)(-?\d+\.\d+)(  )(-?\d+\.\d+)(.*)$/;
    const rGroupsRegex = /\sRGROUPS=\(\d+(\s+\d+)*\)/;

    return this.rawAtomLines.map((line: string, idx: number) => {
      const coordinates = this.coordinates[idx];
      const x = coordinates.x.toFixed(PRECISION);
      const y = coordinates.y.toFixed(PRECISION);

      return line.replace(coordinateRegex, (match, p1, p2, p3, p4, p5) => {
        return p1 + x + p3 + y + p5;
      }).replace(rGroupsRegex, '');
    });
  }
}

