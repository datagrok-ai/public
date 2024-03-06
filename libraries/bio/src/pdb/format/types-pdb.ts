import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {AtomBase, AtomCoordsBase, AtomTerBase} from './types-base';
import {IPdbAtomCoords, IPdbAtomTer} from './types';

export class PdbAtomTer extends AtomTerBase implements IPdbAtomTer {
  constructor(a: AtomBase) {
    super(a, a.number, a.atomElement, a.atomName, a.altLoc, a.resName, a.chain, a.resNumber, a.insCode);
  }
}

/** https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html */
export class PdbAtomCoords extends AtomCoordsBase implements IPdbAtomCoords {
  /**
   * cols: 67-72, SPACE
   * @param segment cols: 73-76
   * @param element cols: 77-78
   * @param charge  cols: 79-80, examples: '1+', '1-'
   */
  constructor(a: AtomCoordsBase,
    public readonly segment: string,
    public readonly element: string,
    public readonly charge: string,
  ) {
    super(a, a.x, a.y, a.z, a.occupancy, a.bFactor);
  }

  static fromStr(src: string): PdbAtomCoords {
    throw new Error('Not implemented');
  }

  override isValid(): boolean {
    return super.isValid() &&
      this.segment.length <= 4 &&
      this.element.length <= 2 &&
      /^()|(\d[+-])$/.test(this.charge);
  }

  override toStr(): string {
    const res = [super.toStr() /* 1-66 */, '      '/* 67-72*/,
      this.segment.padEnd(4) /* 73-76 */,
      this.element.padStart(2) /* 77-78 */,
      this.charge.padEnd(2) /* 79-80 */,
    ].join('');
    if (res.length !== 80)
      throw new Error(`ATOM PDB line of invalid length ${res.length} != 80.`);
    return res;
  }
}
