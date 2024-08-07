// interface for typed arrays, like Float32Array and Uint32Array
export interface ITypedArray {
    length: number;
    [key: number]: any;
}

/** Type for user settings of monomer library set to use. */
export type UserLibSettings = {
  exclude: string[],
  explicit: string[],
}

/** Stores necessary data about atoms of a monomer parsed from Molfile */
export type Atoms = {
    /** Element symbols for monomer's atoms */
    atomTypes: string[],
    /** Cartesian coordiantes of monomer's atoms */
    x: Float32Array,
    y: Float32Array,
    /** V3K atom line may contain keyword args */
    kwargs: string[],
  }

/** Stores necessary data about bonds of a monomer parsed from Molfile */
export type Bonds = {
    /** bond types for all lines of Molfile bond block */
    bondTypes: Uint32Array,
    /** Indices of all atom pairs, indexing starting from 1  */
    atomPairs: number[][],
    /** If a bond has CFG=... keyword argument, it is parsed and sotred as a
     * value of the map, with the key being the bond's index  */
    bondConfiguration: Map<number, number>,
    /** V3K bond line may contain keyword args */
    kwargs: Map<number, string>,
}

/** Metadata associated with the monomer necessary to restore the resulting molfile */
export type MonomerMetadata = {
    /** terminal nodes: 0-th corresponds to the "leftmost" one, 1st, to the "rightmost",
     * e.g. N-terminus and C-terminus in peptides */
    terminalNodes: number[],
    /** r-group nodes: 0-th corresponds to the "leftmost" one, 1st, to the "rightmost" */
    rNodes: number[],
    /** shift from the origin to the next backbone, null for branch monomers */
    backboneShift: number[] | null,
    /** shift from the origin to the next branch, null for branch monomers */
    branchShift: number[] | null
}

export type MolGraph = {
    atoms: Atoms,
    bonds: Bonds,
    meta: MonomerMetadata,
    stereoAtoms?: number[]
}

export type Point = {
    x: number,
    y: number
}

/** Helper structure wrapping common arguments to several functions */
export type LoopVariables = {
    i: number,
    nodeShift: number,
    bondShift: number,
    backbonePositionShift: number[],
    backboneAttachNode: number; // node to which the next backbone is attached
    branchPositionShift: number[],
    branchAttachNode: number,
    flipFactor: number,
    // todo: should we consider representations other than planar?
}

/** Helper structure wrapping common arguments to several functions */
export type LoopConstants = {
    sugar: MolGraph | null,
    phosphate: MolGraph | null,
    seqLength: number,
    atomCount: number,
    bondCount: number,
}

/** Helper structure to simulate pointer to number  */
export type NumberWrapper = {
    value: number | null // null if there is no branch attach node
}
