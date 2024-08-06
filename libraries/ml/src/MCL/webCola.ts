/* eslint-disable max-len */

import Graph from 'graphology';
import forceAtlas2 from 'graphology-layout-forceatlas2';
export function getWebColaLayot(
  is: number[], js: number[], cluster: number[], subCluster: number[]
): {
    embedX: Float32Array;
    embedY: Float32Array;
} {
  const graph = new Graph();
  for (let i = 0; i < cluster.length; i++)
    graph.addNode(cluster[i], {x: Math.random() * 100, y: Math.random() * 100, size: 1});

  for (let i = 0; i < is.length; i++) {
    const p1 = is[i];
    const p2 = js[i];
    const subC1 = subCluster[cluster[p1]];
    const subC2 = subCluster[cluster[p2]];
    const weight = subC1 === subC2 ? 2 : 1;
    graph.addEdge(cluster[is[i]], cluster[js[i]], {weight});
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
    settings: {...forceAtlas2.inferSettings(graph), weighted: true}
  };
  forceAtlas2.inferSettings(graph);
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

  //

  // cluter array contains the indexes of points which are in this big superClueste
  // subCluster contains per index of point, to which subCluster it belongs

  // first, create centroid for each subCluster, and calculate its apparent size
  //   const subClusterSizes = new Map<number, number>();
  //   const interClusterLinks = new Map<number, Map<number, number>>();
  //   let maxInterClusterLinks = 0;
  //   for (let i = 0; i < is.length; i++) {
  //     const c1 = subCluster[cluster[is[i]]];
  //     const c2 = subCluster[cluster[js[i]]];
  //     if (c1 === c2)
  //       continue;
  //     if (!interClusterLinks.has(c1))
  //       interClusterLinks.set(c1, new Map<number, number>());
  //     if (!interClusterLinks.has(c2))
  //       interClusterLinks.set(c2, new Map<number, number>());
  //     const links1 = interClusterLinks.get(c1)!;
  //     const links2 = interClusterLinks.get(c2)!;
  //     links1.set(c2, (links1.get(c2) ?? 0) + 1);
  //     links2.set(c1, (links2.get(c1) ?? 0) + 1);
  //     maxInterClusterLinks = Math.max(maxInterClusterLinks, links1.get(c2)!);
  //     maxInterClusterLinks = Math.max(maxInterClusterLinks, links2.get(c1)!);
  //   }
  //   // calculate the size of each subCluster
  //   for (let i = 0; i < cluster.length; i++) {
  //     const p = cluster[i];
  //     const c = subCluster[p];
  //     if (!subClusterSizes.has(c))
  //       subClusterSizes.set(c, 0);
  //     subClusterSizes.set(c, subClusterSizes.get(c)! + 1);
  //   }

  //   // normilize interClusterLinks
  //   interClusterLinks.forEach((links) => {
  //     links.forEach((count, c) => {
  //       links.set(c, count / maxInterClusterLinks);
  //     });
  //   });


  //   // perform the layout on apperant centroids
  //   // normilaze sizes according to square root of number of points in the cluster
  //   const sqrtNumOfPoints = Math.sqrt(cluster.length);
  //   const normilizedSizes = new Map<number, number>();
  //   subClusterSizes.forEach((size, c) => {
  //     normilizedSizes.set(c, Math.sqrt(size));
  //   });

  //   const maxSize = Array.from(subClusterSizes.values()).reduce((a, b) => Math.max(a, b));

  //   const graph = new Graph();
  //   const clusterNumbers = Array.from(subClusterSizes.keys());
  //   const nodes = clusterNumbers
  //     .map((c, i) => ({index: i, size: normilizedSizes.get(c), x: Math.random() * 10, y: Math.random() * 10}));

  //   for (const node of nodes)
  //     graph.addNode(node.index, {size: node.size, x: node.x, y: node.y});

  //   for (let i = 0; i < clusterNumbers.length; i++) {
  //     if (!interClusterLinks.has(clusterNumbers[i]))
  //       continue;
  //     for (let j = i+1; j < clusterNumbers.length; j++) {
  //       if (!interClusterLinks.get(clusterNumbers[i])!.has(clusterNumbers[j]))
  //         continue;
  //       graph.addEdge(i, j, {weight: interClusterLinks.get(clusterNumbers[i])!.get(clusterNumbers[j])!});
  //     //   links.push({source: i, target: j, weight: interClusterLinks.get(clusterNumbers[i])!.get(clusterNumbers[j])!});
  //     }
  //   }
  //   const scalingRatio = 10;
  //   const settings = {
  //     iterations: 100,
  //     edgeWeightInfluence: 1, // Influence of edge weights on the layout
  //     scalingRatio: scalingRatio, // Scaling ratio to account for node sizes
  //     barnesHutOptimize: true,
  //     barnesHutTheta: 0.5,
  //     adjustSizes: true,
  //     weighted: true
  //   };

  //   forceAtlas2.assign(graph, settings);

  //   //   const layout = new Layout()
  //   //     .size([layoutSize, layoutSize])
  //   //     .nodes(nodes)
  //   //     .links(links)
  //   //     .avoidOverlaps(true)
  //   //     .start(20);

  //   //   for (let i = 0; i < 50; i++) {
  //   //     //@ts-ignore
  //   //     layout.tick();
  //   //   }

  //   const centerX = new Float32Array(nodes.length);
  //   const centerY = new Float32Array(nodes.length);
  //   for (let i = 0; i < nodes.length; i++) {
  //     const node = graph.getNodeAttributes(i);
  //     centerX[i] = node.x;
  //     centerY[i] = node.y;
  //   }
  //   const embedX = new Float32Array(cluster.length);
  //   const embedY = new Float32Array(cluster.length);

  //   // create mapping between subCluster and its index in clusterNumbers
  //   const subClusterIndex = new Map<number, number>();
  //   clusterNumbers.forEach((c, i) => {
  //     subClusterIndex.set(c, i);
  //   });

  //   for (let i = 0; i < cluster.length; i++) {
  //     const c = subCluster[cluster[i]];
  //     const index = subClusterIndex.get(c)!;
  //     const size = normilizedSizes.get(c)!;
  //     embedX[i] = centerX[index] + Math.sin(Math.random() * Math.PI * 2) * size / Math.SQRT2 * Math.random() / maxSize * scalingRatio;
  //     embedY[i] = centerY[index] + Math.cos(Math.random() * Math.PI * 2) * size / Math.SQRT2 * Math.random() / maxSize * scalingRatio;
  //   }

  //   // normilize embeddings
  //   let minX = embedX[0];
  //   let minY = embedY[0];
  //   let maxX = embedX[0];
  //   let maxY = embedY[0];
  //   for (let i = 1; i < cluster.length; i++) {
  //     minX = Math.min(minX, embedX[i]);
  //     minY = Math.min(minY, embedY[i]);
  //     maxX = Math.max(maxX, embedX[i]);
  //     maxY = Math.max(maxY, embedY[i]);
  //   }
  //   let scaleX = maxX - minX;
  //   let scaleY = maxY - minY;
  //   if (scaleX === 0)
  //     scaleX = 1;
  //   if (scaleY === 0)
  //     scaleY = 1;
  //   for (let i = 0; i < cluster.length; i++) {
  //     embedX[i] = (embedX[i] - minX) / scaleX;
  //     embedY[i] = (embedY[i] - minY) / scaleY;
  //   }


//   return {embedX, embedY};
}
