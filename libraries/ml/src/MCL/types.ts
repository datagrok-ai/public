export type MCLOptions = {
    expandFactor: number,
    maxIterations: number,
    inflateFactor: number,
    multFactor: number,
}

export type SparseMatrixObject = {[_: number]: {[_: number]: number}};

export const MCLMethodName = 'MCL';
