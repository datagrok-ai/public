/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { DistanceMatrix } from '@datagrok-libraries/bio/src/trees/distance-matrix';
import { DistanceMetric } from '@datagrok-libraries/bio/src/trees';

export async function getDistanceMatrixForNumerics(values: Float32Array): Promise<Float32Array> {

    return new Promise(function(resolve, reject) {
    const worker = new Worker(new URL('./numerics-distance-worker.ts', import.meta.url));
    worker.postMessage({ values });
    worker.onmessage = ({ data: { error, distanceMatrixData } }): void => {
      error ? reject(error) : resolve(distanceMatrixData);
    };
  });

}

export async function getDistanceMatrixForSequences(sequences: string[]): Promise<Float32Array> {

    return new Promise(function(resolve, reject) {
    const worker = new Worker(new URL('./levenstein-worker.ts', import.meta.url));
    worker.postMessage({ sequences });
    worker.onmessage = ({ data: { error, distanceMatrixData } }): void => {
      error ? reject(error) : resolve(distanceMatrixData);
    };
  });

}
