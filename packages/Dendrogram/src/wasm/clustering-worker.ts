import {getClustersFromDistMatWasm} from './clusterizerWasm';
import {ClusteringWorkerInput} from './clustering-worker-creator';
import {ClusterMatrix} from '@datagrok-libraries/bio/src/trees';


onmessage = (event) => {
  const {distMatArray, n, methodCode}:ClusteringWorkerInput = event.data;
  const data: { error?: any; clusterMatrix?: ClusterMatrix } = {};
  try {
    getClustersFromDistMatWasm(distMatArray, n, methodCode)
      .then((clusterMatrix) => {
        data.clusterMatrix = clusterMatrix;
        postMessage(data);
      });
  } catch (e) {
    data.error = e;
    postMessage(data);
  }
};
