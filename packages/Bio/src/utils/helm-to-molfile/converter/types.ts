export type Bond = {
  /** Global (for complex polymer) or local (for simple polymer) monomer index, starting from 0 */
  monomerIdx: number,
  /** RGroup id, starting from 1  */
  rGroupId: number
}

/** Position of a node in the connection list / bond block  */
export type PositionInBonds = {
  bondLineIdx: number,
  nodeIdx: number,
}

export type MonomerMap = { position: number, symbol: string, atoms: number[], bonds: number[] };

export class MolfileWithMap {
  constructor(
    public readonly molfile: string,
    public readonly monomers: MonomerMap[]
  ) {}

  static empty() { return new MolfileWithMap('', []); }
}
