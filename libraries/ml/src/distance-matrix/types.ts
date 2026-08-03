export enum DistanceAggregationMethods {
    EUCLIDEAN = 'EUCLIDEAN',
    MANHATTAN = 'MANHATTAN',
  };

export type DistanceAggregationMethod = keyof typeof DistanceAggregationMethods;

/** Per-cluster representativeness. All three arrays are aligned with the cluster's input member
 * list (the row indices passed in `clusters`). `meanDistances[k]` is member k's average distance to
 * the other members (0 for a singleton); `ranks[k]` is its 1-based rank by ascending mean distance
 * (rank 1 = medoid), ties broken by ascending member index. */
export type ClusterRepresentatives = {
  members: number[];
  meanDistances: number[];
  ranks: number[];
};
