/* eslint-disable max-len */

import Graph from 'graphology';
import forceAtlas2 from 'graphology-layout-forceatlas2';

// note: subCluster is of size nRows and contains for each point which subcluster it belongs to
// note: cluster is of size of this concrete supercluster, and contains the indexes of the points in that supercluster
export function getWebColaLayot(
  cluster: number[], clusterConnections: {i: number[], j: number[], v: number[]}, _subCluster: number[]
): {
    embedX: Float32Array;
    embedY: Float32Array;
} {
  const graph = new Graph();
  for (let i = 0; i < cluster.length; i++)
    graph.addNode(cluster[i], {x: Math.random() * 20000, y: Math.random() * 20000});

  for (let it = 0; it < clusterConnections.i.length; it++) {
    const i = clusterConnections.i[it];
    const j = clusterConnections.j[it];
    const v = clusterConnections.v[it];
    if (i === j || v <= 0)
      continue;
    const weight = _subCluster[i] === _subCluster[j] ? 2 : 1;
    graph.addEdge(i, j, {weight: v * weight});
  }

  const settings = {
    iterations: 100,
    // edgeWeightInfluence: 1, // Influence of edge weights on the layout
    // scalingRatio: 1, // Scaling ratio to account for node sizes
    // //barnesHutOptimize: true,
    // //barnesHutTheta: 0.5,
    // adjustSizes: true,
    // weighted: true,
    // strongGravityMode: true,
    getEdgeWeight: 'weight',
    settings: {...forceAtlas2.inferSettings(graph), weighted: true, edgeWeightInfluence: 1}
  };
  forceAtlas2.assign(graph, settings);

  const embedX1 = new Float32Array(cluster.length);
  const embedY1 = new Float32Array(cluster.length);
  for (let i = 0; i < cluster.length; i++) {
    const node = graph.getNodeAttributes(cluster[i]);
    embedX1[i] = node.x;
    embedY1[i] = node.y;
  }

  let minX = embedX1[0];
  let minY = embedY1[0];
  let maxX = embedX1[0];
  let maxY = embedY1[0];
  for (let i = 1; i < cluster.length; i++) {
    minX = Math.min(minX, embedX1[i]);
    minY = Math.min(minY, embedY1[i]);
    maxX = Math.max(maxX, embedX1[i]);
    maxY = Math.max(maxY, embedY1[i]);
  }
  let scaleX = maxX - minX;
  let scaleY = maxY - minY;
  if (scaleX === 0)
    scaleX = 1;
  if (scaleY === 0)
    scaleY = 1;
  for (let i = 0; i < cluster.length; i++) {
    embedX1[i] = (embedX1[i] - minX) / scaleX;
    embedY1[i] = (embedY1[i] - minY) / scaleY;
  }

  return {embedX: embedX1, embedY: embedY1};
}
