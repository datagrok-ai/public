import {MolfileHandlerBase} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler-base';
import {MolfileAtoms} from './mol-atoms';

export class MolfileAtomsV2K extends MolfileAtoms {
  constructor(molfileHandler: MolfileHandlerBase) {
    super();
    this.rawAtomLines = molfileHandler.getAtomLines();
    this.coordinates = this.rawAtomLines.map((line: string) => {
      const x = parseFloat(line.substring(0, 10));
      const y = parseFloat(line.substring(10, 20));
      return {x, y};
    });
  }

  get atomLines(): string[] {
    return this.rawAtomLines.map((line: string, idx: number) => {
      const coordinates = this.coordinates[idx];
      const x = coordinates.x.toFixed(4).padStart(10, ' ');
      const y = coordinates.y.toFixed(4).padStart(10, ' ');
      return `${x}${y}${line.substring(20)}`;
    });
  }
}

