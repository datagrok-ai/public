import { DistanceMatrix } from "@datagrok-libraries/bio/src/trees/distance-matrix";

onmessage = (event) => {
  const { values } = event.data;

  const data: { error?: any; distanceMatrixData?: Float32Array } = {};
  try {
    const distanceMatrix = DistanceMatrix.calc(values, (a: number, b: number) => Math.abs(a - b));
    data.distanceMatrixData = distanceMatrix.data;
  } catch (e) {
    data.error = e;
  }
  postMessage(data);
};
