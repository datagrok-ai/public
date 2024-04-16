export type SparseMatrix = {
    i: Uint32Array;
    j: Uint32Array;
    distance: Float32Array;
}

export type Matrix = Array<Uint32Array> | Array<Int32Array> | Array<Float32Array>;

export enum MatrixMatrixOpType {
    ADD='ADD',
    SUB='SUB',
    MULT='MULT',
}
export enum MatrixOpType {
    SQUARE='SQUARE',
    INVERSE='INVERSE',
    TRANSPOSE='TRANSPOSE',
    NORM = 'NORM',
    COLUMN_NORM = 'COLUMN_NORM',
}

export enum MatrixScalarOpType {
    SCALARMULT='SCALARMULT',
    SCALARADD='SCALARADD',
    SCALARPOW='SCALARPOW',
}
