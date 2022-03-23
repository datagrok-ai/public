import {Options, Coordinates, Vectors, DistanceMetric} from '@datagrok-libraries/utils/src/type-declarations';
import {
  calculateEuclideanDistance,
  calcDistanceMatrix,
  fillRandomMatrix,
  vectorAdd,
} from '@datagrok-libraries/utils/src/vector-operations';
import {randomInt} from '@datagrok-libraries/utils/src/random';

/**
 * Implements stochastic proximity embedding.
 *
 * @export
 * @class SPEBase
 * @link doi:10.1016/S1093-3263(03)00155-4
 */
export class SPEBase {
  protected static dimension = 2;
  protected steps: number;
  protected cycles: number;
  protected cutoff: number;
  protected lambda: number;
  protected dlambda: number;
  protected lambda2: number;
  protected dlambda2: number;
  protected epsilon: number;
  protected distanceFunction: DistanceMetric;
  protected distance: Coordinates;

  /**
   * Creates an instance of SPEBase.
   * @param {Options} [options] Options to pass to the constructor.
   * @memberof SPEBase
   */
  constructor(options?: Options) {
    this.steps = options?.steps ?? 0;
    this.cycles = options?.cycles ?? 1e6;
    // Select a cutoff distance {cutoff} and ...
    this.cutoff = options?.cutoff ?? 0;
    // ... an initial learning rate {lambda} > 0
    this.lambda = options?.lambda ?? 2.0;
    this.dlambda = options?.dlambda ?? 0.01;
    this.lambda2 = this.lambda/2.;
    this.dlambda2 = this.dlambda/2.;
    this.epsilon = options?.epsilon ?? 1e-10;
    // eslint-disable-next-line brace-style
    this.distanceFunction = options?.distance ?? calculateEuclideanDistance;
    this.distance = [];
  }

  /**
   * Initializes distance matrix.
   *
   * @protected
   * @param {Vectors} vectors Input vectors to calculate distance between.
   * @memberof SPEBase
   */
  protected initDistance(vectors: Vectors) {
    this.distance = calcDistanceMatrix(vectors, this.distanceFunction);
  }

  /**
   * Calculates distance between the two vectors given.
   *
   * @param {Vectors} vectors Set of vectors to calculate distances between.
   * @param {number} index1 Index of the first vector of the pair.
   * @param {number} index2 Index of the second vector of the pair.
   * @return {number} Distance between these two vectors.
   */
  protected calcDistance(vectors: Vectors, index1: number, index2: number): number {
    return this.distance[index1][index2];
  }

  /**
   * Embeds the vectors given into a two-dimensional space.
   *
   * @param {Vectors} vectors D-dimensional coordinates.
   * @return {Coordinates} SPE coordinates in D space.
   */
  public embed(vectors: Vectors): Coordinates {
    const nItems = vectors.length;
    const areaWidth = 40;
    // Initialize the D-dimensional coordinates of the N points.
    const coordinates = fillRandomMatrix(nItems, SPEBase.dimension, areaWidth);

    let lambda2 = this.lambda2;

    if (this.steps == 0) {
      this.steps = vectors.length-1;
    }

    this.initDistance(vectors);

    for (let cycle = 0; cycle < this.cycles; ++cycle) {
      for (let step = 0; step < this.steps; ++step) {
        // Select two points, i and j, at random, ...
        const i = randomInt(nItems); let j = randomInt(nItems);
        while (i == j) j = randomInt(nItems);

        const rowi = coordinates[i]; const rowj = coordinates[j];

        // ... retrieve (or evaluate) their proximity in the input space, rij and ...
        const r = this.calcDistance(vectors, i, j);
        // ... compute their Euclidean distance on the D-dimensional map, dij.
        const d = calculateEuclideanDistance(rowi, rowj);

        // If rij <= rc, or if rij > rc and dij < rij ...
        if ((this.cutoff == 0) || (r <= this.cutoff) || (d < r)) {
          const multiplier = lambda2*(r-d)/(d+this.epsilon);
          // ... update the coordinates xi and xj.
          const diffIJ = vectorAdd(rowi, rowj, -1);
          coordinates[i] = vectorAdd(rowi, diffIJ, multiplier);
          coordinates[j] = vectorAdd(rowj, diffIJ, -multiplier);
        }
      }
      // Decrease the learning rate {lambda} by a prescribed {dlambda}.
      lambda2 -= this.dlambda2;
      if (lambda2 <= 0.) {
        break;
      }
    }
    return coordinates;
  }
}

/**
 * Implements modified stochastic proximity embedding.
 *
 * @export
 * @class PSPEBase
 * @link doi:10.1016/S1093-3263(03)00155-4
 */
export class PSPEBase extends SPEBase {
  /**
   * Embeds the vectors given into a two-dimensional space using a modified update rule.
   *
   * @param {Vectors} vectors D-dimensional coordinates.
   * @return {Coordinates} SPE coordinates in D space.
   */
  public embed(vectors: Vectors): Coordinates {
    const nItems = vectors.length;
    const areaWidth = 40;
    //  Initialize the D-dimensional coordinates of the N points.
    const coordinates = fillRandomMatrix(nItems, PSPEBase.dimension, areaWidth);
    let lambda = this.lambda;

    this.initDistance(vectors);

    for (let cycle = 0; cycle < this.cycles; ++cycle) {
      // Select a point, i, at random (pivot).
      const i: number = randomInt(nItems);
      const rowi = coordinates[i];

      // For every point j != i ...
      for (let j = 0; j < nItems; ++j) {
        if (i == j) continue;
        const rowj = coordinates[j];
        // ... retrieve (or evaluate) its proximity to i in the input space, rij ...
        const r = this.calcDistance(vectors, i, j);
        // ... and compute their Euclidean distance on the D-dimensional map, dij.
        const d = calculateEuclideanDistance(rowi, rowj);
        // If rij <= rc, or if rij > rc and dij < rij ...
        if ((this.cutoff == 0) || (r <= this.cutoff) || (d < r)) {
          const multiplier = lambda*(r-d)/(d+this.epsilon);
          const diffIJ = vectorAdd(rowi, rowj, -1);
          // ... update the coordinates xj.
          coordinates[j] = vectorAdd(rowj, diffIJ, -multiplier);
        }
      }
      // Decrease the learning rate {lambda} by a prescribed {dlambda}.
      lambda -= this.dlambda;
      if (lambda <= 0.) {
        break;
      }
    }
    return coordinates;
  }
}

/**
 * Implements modified stochastic proximity embedding.
 *
 * @export
 * @class OriginalSPE
 * @link doi:10.1002/jcc.10234
 */
export class OriginalSPE extends SPEBase {
  /**
   * Embeds the vectors given into a two-dimensional space using a modified update rule.
   *
   * @param {Vectors} vectors D-dimensional coordinates.
   * @return {Coordinates} SPE coordinates in D space.
   */
  protected radiusPercent: number;
  protected maxDistance: number;
  protected maxDistanceSteps: number;

  constructor(options?: Options) {
    super(options);
    this.cycles = options?.cycles ?? 1e3;
    this.steps = options?.steps ?? 100000;
    this.radiusPercent = options?.radiusPercent ?? 1.0;
    this.maxDistance = options?.maxDistance ?? null;
    this.maxDistanceSteps = options?.maxDistanceSteps ?? null;
  }

  public embed(vectors: Vectors): Coordinates {
    const nItems = vectors.length;
    const areaWidth = 40;
    // Initialize the D-dimensional coordinates of the N points.
    const coordinates = fillRandomMatrix(nItems, OriginalSPE.dimension, areaWidth);

    this.initDistance(vectors);

    if (this.maxDistanceSteps == null) {
      this.maxDistanceSteps = nItems * Math.floor((nItems - 1) / 2);
    }
    if (this.maxDistance == null) {
      this.maxDistance = -1e37;
      for (let n = 0; n < this.maxDistanceSteps; n++) {
        const i = randomInt(nItems); let j = randomInt(nItems);
        while (i == j) j = randomInt(nItems);

        const d = this.calcDistance(vectors, i, j);
        if (d > this.maxDistance) {
          this.maxDistance = d;
        }
      }
    }

    let lambda = this.lambda;
    const radius = (this.radiusPercent == 0.0) ? this.maxDistance : this.maxDistance * this.radiusPercent;

    for (let cycle = 0; cycle < this.cycles; ++cycle) {
      for (let step = 0; step < this.steps; ++step) {
        // Select two points, i and j, at random, ...
        const i = randomInt(nItems); let j = randomInt(nItems);
        while (i == j) j = randomInt(nItems);

        const rowi = coordinates[i]; const rowj = coordinates[j];

        // ... retrieve (or evaluate) their proximity in the input space, rij and ...
        const r = this.calcDistance(vectors, i, j);
        // ... compute their Euclidean distance on the D-dimensional map, dij.
        const d = calculateEuclideanDistance(rowi, rowj);

        if ((r <= radius) || (d < r)) {
          const multiplier = lambda * 0.5 * (r - d) / (d + this.epsilon);
          // ... update the coordinates xi and xj.
          const diffIJ = vectorAdd(rowi, rowj, -1);
          coordinates[i] = vectorAdd(rowi, diffIJ, multiplier);
          coordinates[j] = vectorAdd(rowj, diffIJ, -multiplier);
        }
      }
      lambda -= ((this.lambda - this.dlambda) / (this.cycles - 1.0)); ;
      if (lambda < this.dlambda) {
        break;
      }
    }
    return coordinates;
  }
}
