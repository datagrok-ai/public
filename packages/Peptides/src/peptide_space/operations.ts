import {Coordinates} from './declarations';

export function transpose(matrix: Coordinates): Coordinates {
  return matrix[0].map((col, i) => matrix.map((row) => row[i]));
}

export function substractVectors(v1: number[], v2: number[]): number[] {
  return v1.map((v, i) => v-v2[i]);
}

export function addVectors(v1: number[], v2: number[]): number[] {
  return v1.map((v, i) => v+v2[i]);
}

export function multiplyVector(vec: number[], n: number): number[] {
  return vec.map((v, i) => v*n);
}

export function calculateEuclideanDistance(p: number[], q: number[]): number {
  if (p.length != q.length) throw new Error('Arrays have different dimensions.');
  const subtracted = substractVectors(p, q);
  const powered = subtracted.map((e) => Math.pow(e, 2));
  const sum = powered.reduce((total, current) => total + current, 0);
  return Math.sqrt(sum);
};
