/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { DistanceMatrix } from '@datagrok-libraries/bio/src/trees/distance-matrix';

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

export async function getCombinedDistanceMatrix(columns: DG.Column[], method: 'euclidean' | 'manhattan' = 'euclidean') {
    const promises: Promise<Float32Array>[] = [];
    
    columns.forEach(col => {
        if(col.type === DG.TYPE.FLOAT || col.type === DG.TYPE.INT) {
            promises.push(getDistanceMatrixForNumerics(col.toList() as unknown as Float32Array));
        }
        else if(col.semType ==='Macromolecule') {
            promises.push(getDistanceMatrixForSequences(col.toList() as unknown as string[]));
        }
    });

    // return the array of resolved promises values after all are resolved
    return Promise.all(promises).then((values) => {
        if (values.length === 1)
            return new DistanceMatrix(values[0]);
        else
            return DistanceMatrix.combineDistanceMatrices(values.map(v => new DistanceMatrix(v)), method);
    }).catch((error) => {
        grok.shell.error(error);
    });
}