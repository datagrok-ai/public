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

/** Cap group information for an R-group */
export type CapGroupInfo = {
  /** Extracted element string (e.g. 'H', 'O', 'C=C') */
  element: string,
  /** Raw cap group SMILES (e.g. '[H][*:1]', 'O[*:2]', 'C=C[*:3]') */
  smiles: string,
  /** Whether the cap is a single atom (valid element symbol) */
  isSimple: boolean,
  /** Number of R group, to handle cases where its not sorted */
  rGroupId: number,
}

