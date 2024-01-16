export enum DistanceAggregationMethods {
    EUCLIDEAN = 'EUCLIDEAN',
    MANHATTAN = 'MANHATTAN',
  };

export type DistanceAggregationMethod = keyof typeof DistanceAggregationMethods;
