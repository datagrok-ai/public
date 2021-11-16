import * as vec from 'vectorious';

export import Matrix = vec.NDArray;

export type Options = {[name: string]: any};
export type DistanceMetric = (v1: any, v2: any) => (number);
export type Coordinates = Array<Float32Array>;
export type Vectors = Array<any>;
