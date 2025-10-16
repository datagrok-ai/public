export type NumericArray = Float32Array | Float64Array | Int32Array | Uint32Array;
export type Sense = "min" | "max";
export type ParetoLabel = "optimal" | "non-optimal";

/**
 * Optimized Pareto mask computation (coordinate-wise input)
 * @param data Array of numeric arrays where each array contains one dimension for all points
 * @param sense Array of "min" or "max" for each dimension
 */
export function paretoMaskFromCoordinates(data: NumericArray[], sense: Sense[], nPoints: number): ParetoLabel[] {
  const nDims = data.length;
  if (nDims === 0) return [];
  //const nPoints = data[0].length;

  if (sense.length !== nDims) {
    throw new Error("Sense array length must match number of dimensions");
  }

  // Transpose data: convert from coordinate-wise to points
  const points: NumericArray[] = Array(nPoints)
    .fill(0)
    .map((_, i) => {
      const point = new Float64Array(nDims);
      for (let d = 0; d < nDims; d++) {
        point[d] = data[d][i];
      }
      return point;
    });

  // Reuse optimized Pareto mask function
  return paretoMaskOptimized(points, sense, nPoints);
}

/**
 * Optimized Pareto mask function for points
 */
function paretoMaskOptimized(data: NumericArray[], sense: Sense[], nPoints: number): ParetoLabel[] {
  //const nPoints = data.length;
  if (nPoints === 0) return [];
  const nDims = data[0].length;

  const points = data.map((p, i) => ({ point: p, index: i }));

  // Sort by first dimension according to sense
  points.sort((a, b) => (sense[0] === "min" ? a.point[0] - b.point[0] : b.point[0] - a.point[0]));

  const mask: ParetoLabel[] = Array(nPoints).fill("optimal");
  const paretoFront: NumericArray[] = [];

  for (const { point, index } of points) {
    let dominated = false;

    for (const frontPoint of paretoFront) {
      let dominates = true;
      let strictlyBetter = false;

      for (let d = 0; d < nDims; d++) {
        const a = frontPoint[d];
        const b = point[d];
        const s = sense[d];

        if (s === "min") {
          if (a > b) dominates = false;
          if (a < b) strictlyBetter = true;
        } else {
          if (a < b) dominates = false;
          if (a > b) strictlyBetter = true;
        }
      }

      if (dominates && strictlyBetter) {
        dominated = true;
        break;
      }
    }

    if (dominated) {
      mask[index] = "non-optimal";
    } else {
      paretoFront.push(point);
    }
  }

  return mask;
}

const x = new Float32Array([0, 0, 0, 1, 2, 2, 4, 4.5]);
const y = new Float32Array([4, 2, 0, 1, 0, 3, 1, 0]);

// Example usage (coordinate-wise input):
const data: NumericArray[] = [
  x, // first dimension
  y, // second dimension
];

const sense: Sense[] = ["max", "max"];

const mask = paretoMaskFromCoordinates(data, sense, x.length);

for (let i = 0; i < mask.length; ++i)
    console.log('(', x[i], ',', y[i], ') <--> ', mask[i]);

