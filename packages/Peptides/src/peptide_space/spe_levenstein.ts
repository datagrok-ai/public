//import {Memoize,MemoizeExpiring} from 'typescript-memoize';
import * as fl from 'fastest-levenshtein';
import {SPEBase, Vectors, Coordinates} from './spe';
import {CachedCalculator} from './cache';

export class SequenceSPE extends SPEBase {
  /**
     * distanceFunction
     * @param {any} item1 First vector of the pair.
     * @param {any} item2 Second vector of the pair.
     * @return {number} Distance between these two vectors.
     */
  protected distanceFunction(item1: any, item2: any): number {
    return fl.distance(item1, item2);
  }
}

class SPECachedCalculator extends CachedCalculator {
  protected makeKey(vectors: Vectors, index1: number, index2: number): string {
    return index1.toString()+';'+index2.toString();
  }
}

export class SequenceSPECached extends SequenceSPE {
  protected cache: SPECachedCalculator;

  constructor(
    steps: number,
    cycles: number,
    cutoff: number,
    lambda: number,
    dlambda: number,
    epsilon: number = 1e-10,
  ) {
    super(steps, cycles, cutoff, lambda, dlambda, epsilon);
    this.cache = new SPECachedCalculator();
  }

  protected calcDistanceDirect(vectors: Vectors, index1: number, index2: number): number {
    return super.distanceFunction(vectors[index1], vectors[index2]);
  }

  /**
   * calcDistance
   * @param {Vectors} vectors Set of vectors to calculate distances between.
   * @param {number} index1 Index of the first vector of the pair.
   * @param {number} index2 Index of the second vector of the pair.
   * @return {number} Distance between these two vectors using cached calls.
   */
  protected calcDistance(vectors: Vectors, index1: number, index2: number): number {
    return this.cache.calcCached(this.calcDistanceDirect, vectors, index1, index2);
  }

  /**
   * embed
   * @param {Vectors} vectors D-dimensional coordinates.
   * @return {Coordinates} SPE coordinates in D space using cached calls.
   */
  public embed(vectors: Vectors): Coordinates {
    this.cache.resetCache();
    return super.embed(vectors);
  }
}
