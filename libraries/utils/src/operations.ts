import * as vec from 'vectorious';
import {Matrix, Coordinates} from './type_declarations';

export function transpose(matrix: Matrix): Coordinates {
  return matrix.transpose().toArray();
}

export function calculateEuclideanDistance(p: vec.NDArray, q: vec.NDArray): number {
  return Math.sqrt(vec.subtract(p, q).pow(2).sum());
}

export function fillRandomMatrix(dimension1: number, dimension2: number, scale: number = 1.): Matrix {
  return vec.random(dimension1, dimension2).scale(scale);
}