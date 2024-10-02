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

