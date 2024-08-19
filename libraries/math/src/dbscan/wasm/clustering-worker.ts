import {dbscan} from './dbscan';

onmessage = (event) => {
  const {embedX, embedY, epsilon, minPts} = event.data;
  const data: { error?: any; clusters?: Int32Array } = {};
  (async () => {
    try {
      const clusters = await dbscan(embedX, embedY, epsilon, minPts)
      data.clusters = clusters;
      postMessage(data);
    } catch (e) {
      data.error = e;
      postMessage(data);
    }
  })();
};
