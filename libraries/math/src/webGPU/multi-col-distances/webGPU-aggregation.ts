

function euclideanAggregationWgsl(arraySize: number) {
  return `
        var sum = 0.0;
        for (var i = 0u; i < ${arraySize}; i = i + 1u) {
            sum = sum + distances[i] * distances[i] * computeInfo.weights[i] * computeInfo.weights[i];
        }
        return sqrt(sum);
    `;
};

function manhattanAggregationWgsl(arraySize: number) {
  return `
        var sum = 0.0;
        for (var i = 0u; i < ${arraySize}; i = i + 1u) {
            sum = sum + abs(distances[i]) * computeInfo.weights[i];
        }
        return sum;
    `;
}


export enum WEBGSLAGGREGATION {
    EUCLIDEAN = 'EUCLIDEAN',
    MANHATTAN = 'MANHATTAN'
}

export const WEBGSLAGGREGATIONFUNCTIONS: {[key in WEBGSLAGGREGATION]: (arraySize: number) => string} = {
  [WEBGSLAGGREGATION.EUCLIDEAN]: euclideanAggregationWgsl,
  [WEBGSLAGGREGATION.MANHATTAN]: manhattanAggregationWgsl
};
