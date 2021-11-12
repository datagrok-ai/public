export type Coordinates = Array<Array<number>>;
export type Vectors = Array<any>;

function substractVectors(v1: number[], v2: number[]): number[] {
  return v1.map((v, i) => v-v2[i]);
}

function addVectors(v1: number[], v2: number[]): number[] {
  return v1.map((v, i) => v+v2[i]);
}

function multiplyVector(vec: number[], n: number): number[] {
  return vec.map((v, i) => v*n);
}

const calculateEuclideanDistance = (p: number[], q: number[]) => {
  if (p.length != q.length) throw new Error('Arrays have different dimensions.');
  const subtracted = substractVectors(p, q);
  const powered = subtracted.map((e) => Math.pow(e, 2));
  const sum = powered.reduce((total, current) => total + current, 0);
  return Math.sqrt(sum);
};

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
  protected distanceFunction: Function;

  constructor(options: {[name: string]: any} = {
    steps: 1000,
    cycles: 1000,
    cutoff: 0,
    lambda: 2.0,
    dlambda: 0.01,
    epsilon: 1e-10,
    distance: Function,
  }) {
    this.steps = options.steps;
    this.cycles = options.cycles;
    // Select a cutoff distance {cutoff} and ...
    this.cutoff = options.cutoff;
    // ... an initial learning rate {lambda} > 0
    this.lambda = options.lambda;
    this.dlambda = options.dlambda;
    this.lambda2 = options.lambda/2.;
    this.dlambda2 = options.dlambda/2.;
    this.epsilon = options.epsilon;
    this.distanceFunction = options.distance;
  }

  /**
   * calcDistance
   * @param {Vectors} vectors Set of vectors to calculate distances between.
   * @param {number} index1 Index of the first vector of the pair.
   * @param {number} index2 Index of the second vector of the pair.
   * @return {number} Distance between these two vectors.
   */
  protected calcDistance(vectors: Vectors, index1: number, index2: number): number {
    return this.distanceFunction(vectors[index1], vectors[index2]);
  }

  /**
   * embed
   * @param {Vectors} vectors D-dimensional coordinates.
   * @return {Coordinates} SPE coordinates in D space.
   */
  public embed(vectors: Vectors): Coordinates {
    const nItems = vectors.length;
    const areaWidth = 40;
    // Initialize the D-dimensional coordinates of the N points.
    const coordinates: Coordinates = new Array(nItems).fill(0).map(
      () => Array.from({length: SPEBase.dimension}, () => Math.floor(Math.random() * areaWidth)),
    );

    const _randomInt = (range: number) => (Math.floor(Math.random() * range));
    let lambda2 = this.lambda2;

    for (let cycle = 0; cycle < this.cycles; ++cycle) {
      for (let step = 0; step < this.steps; ++step) {
        // Select two points, i and j, at random, ...
        const i = _randomInt(nItems); let j = _randomInt(nItems);
        while (i == j) j = _randomInt(nItems);

        // ... retrieve (or evaluate) their proximity in the input space, rij and ...
        const r = this.calcDistance(vectors, i, j);
        // ... compute their Euclidean distance on the D-dimensional map, dij.
        const d = calculateEuclideanDistance(coordinates[i], coordinates[j]);

        // If rij <= rc, or if rij > rc and dij < rij ...
        if ((this.cutoff == 0) || (r <= this.cutoff) || (d < r)) {
          const multiplier = lambda2*(r-d)/(d+this.epsilon);
          // ... update the coordinates xi and xj.
          coordinates[i] = addVectors(
            coordinates[i],
            multiplyVector(substractVectors(coordinates[i], coordinates[j]), multiplier),
          );
          coordinates[j] = addVectors(
            coordinates[j],
            multiplyVector(substractVectors(coordinates[j], coordinates[i]), multiplier),
          );
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

export class PSPEBase extends SPEBase {
  /**
   * embed
   * @param {Vectors} vectors D-dimensional coordinates.
   * @return {Coordinates} SPE coordinates in D space.
   */
  public embed(vectors: Vectors): Coordinates {
    const nItems = vectors.length;
    const areaWidth = 40;
    //  Initialize the D-dimensional coordinates of the N points.
    const coordinates: Coordinates = new Array(nItems).fill(0).map(
      () => Array.from({length: PSPEBase.dimension}, () => Math.floor(Math.random() * areaWidth)),
    );

    const _randomInt = (range: number) => (Math.floor(Math.random() * range));
    let lambda = this.lambda;

    for (let cycle = 0; cycle < this.cycles; ++cycle) {
      // Select a point, i, at random (pivot).
      const i = _randomInt(nItems);

      // For every point j != i ...
      for (let j = 0; j < nItems; ++j) {
        if (i == j) continue;
        // ... retrieve (or evaluate) its proximity to i in the input space, rij ...
        const r = this.calcDistance(vectors, i, j);
        // ... and compute their Euclidean distance on the D-dimensional map, dij.
        const d = calculateEuclideanDistance(coordinates[i], coordinates[j]);
        // If rij <= rc, or if rij > rc and dij < rij ...
        if ((this.cutoff == 0) || (r <= this.cutoff) || (d < r)) {
          const multiplier = lambda*(r-d)/(d+this.epsilon);
          // ... update the coordinates xj.
          coordinates[j] = addVectors(
            coordinates[j],
            multiplyVector(substractVectors(coordinates[j], coordinates[i]), multiplier),
          );
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
