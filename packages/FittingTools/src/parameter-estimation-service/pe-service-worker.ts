import {PeServiceWorkerBase} from "./pe-service-worker-base";
import {fit, FitErrorModel, sigmoid} from "@datagrok-libraries/statistics/src/parameter-estimation/fit-curve";
export class PeServiceWorker extends PeServiceWorkerBase {
  fit(x: number[], y: number[], counts: number[]): number [] {
    const paramCount = 4;
/*
    if (counts.length > 0) {
      console.log(x);
      console.log(y);
      console.log(counts);
      console.log('Total Length ' + x.length);
    }
*/

    const ar = new Array<number>(paramCount*counts.length);
    let offset = 0;
    let offsetParam = 0;
    for (let n = 0; n < counts.length; ++n) {
      //console.log('offset: ' + offset + " count " + counts[n] + " x.length = " + x.length);
      const segX = x.slice(offset, offset + counts[n]);
      const segY = y.slice(offset, offset + counts[n]);

      const segYCopy = segY.slice(0);
      segYCopy.sort((a, b) => a - b);
      //console.log(segYCopy);
      const minY = segYCopy[0];
      const maxY = segYCopy[segYCopy.length -1];
      const idxMed = Math.ceil(segYCopy.length / 2);
      const medY = segYCopy[idxMed];
      let xAtMedY = -1;
      for (let k = 0; k < segY.length; ++k) {
        if (segY[k] === medY) {
          xAtMedY = segX[k];
          break;
        }
      }

      //console.log(segX);
      //console.log(segY);
      offset += counts[n];
      let fitRes = fit({x: segX, y: segY}, [maxY, 1.2, xAtMedY, minY], sigmoid, FitErrorModel.Constant);
      let params = fitRes.parameters;
      //console.log('param: ' + params);
      offsetParam = paramCount * n;
      for (let i = 0; i < params.length; ++i) {
        ar[offsetParam + i] = params[i];
      }
    }
    //console.log('Finished Worker ' + ar + " " + counts.length);
    return ar;
  }
}
