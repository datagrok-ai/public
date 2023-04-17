import { DistanceMatrix } from "@datagrok-libraries/bio/src/trees/distance-matrix";

onmessage = (event) => {
  const { values } = event.data;

  const distanceFunc = (a: number, b: number) => {
    return Math.abs(a - b);
  };

  const data: { error?: any; distanceMatrix?: DistanceMatrix } = {};
  try {
    const distanceMatrix = DistanceMatrix.calc(values, distanceFunc);
    data.distanceMatrix = distanceMatrix;
  } catch (e) {
    data.error = e;
  }
  postMessage(data);
};
