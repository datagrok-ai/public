import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {AtomBase, AtomCoordsBase, AtomTerBase} from './types-base';
import {IPdbAtomTer, IPdbqtAtomCoords, IPdbqtAtomTer} from './types';
import {PdbAtomCoords, PdbAtomTer} from './types-pdb';

export type VinaResultType = {
  affinity: number,
  lbDistFromBest: number,
  ubDistFromBest: number
};

export class PdbqtAtomTer extends AtomTerBase implements IPdbqtAtomTer {
  constructor(a: AtomBase) {
    super(a, a.number, a.atomElement, a.atomName, a.altLoc, a.resName, a.chain, a.resNumber, a.insCode);
  }

  public toPdb(): IPdbAtomTer {
    return new PdbAtomTer(this);
  }

  static fromStr(line: string): PdbqtAtomTer {
    return new PdbqtAtomTer(AtomBase.fromStr(line));
  }
}

/** https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html */
export class PdbqtAtomCoords extends AtomCoordsBase implements IPdbqtAtomCoords {
  /**
   * cols: 67-70, SPACE
   * @param partCharge cols: 71-76 f:+1.3
   * @param element    cols: 78-79, right
   */
  constructor(a: AtomCoordsBase,
    public readonly partCharge: number,
    public readonly element: string,
  ) {
    super(a, a.x, a.y, a.z, a.occupancy, a.bFactor);
  }

  static fromStr(line: string): PdbqtAtomCoords {
    return new PdbqtAtomCoords(
      AtomCoordsBase.fromStr(line),
      /* partCharge */ parseFloat(line.slice(71 - 1, 76)),
      `${line.slice(78 - 1, 79).trim()}`,
    );
  }

  /** cols: 71-76 f:+5.3 */ get partChargeStr(): string {
    return `${this.partCharge >= 0 ? '+' : ''}${this.partCharge.toFixed(3)}`;
  }

  // /** cols: 73-76, left */ segment: string;
  isValid(): boolean {
    return super.isValid() &&
      this.partChargeStr.length <= 6 &&
      this.element.length <= 2;
  }

  toPdb(): PdbAtomCoords {
    const res = new PdbAtomCoords(this, '', this.element, '');
    return res;
  }

  toStr(): string {
    const res = [super.toStr() /* 1-66 */, '    ' /* 67-70 */,
      this.partChargeStr.padStart(6) /* 71-76*/, ' ' /* 77 */,
      this.element.padEnd(2) /* 78-79 */,
    ].join('');
    if (res.length !== 79)
      throw new Error(`ATOM Pdbqt line of invalid length ${res.length} != 79.`);
    return res;
  }
}
