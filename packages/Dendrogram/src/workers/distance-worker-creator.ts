/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { DistanceMatrix } from '@datagrok-libraries/bio/src/trees/distance-matrix';

export async function getDistanceMatrixForNumerics(data: Float32Array) {
    
    return new Promise(function(resolve, reject) {
    const worker = new Worker(new URL('./numerics-distance-worker.ts', import.meta.url));
    worker.postMessage({data});
    worker.onmessage = ({data: {error, distanceMatrix}}): void => {
      error ? reject(error) : resolve(distanceMatrix);
    };
  });

}

export async function getDistanceMatrixForSequences(data: string[]) {
    
    return new Promise(function(resolve, reject) {
    const worker = new Worker(new URL('./levenstein-worker.ts', import.meta.url));
    worker.postMessage({data});
    worker.onmessage = ({data: {error, distanceMatrix}}): void => {
      error ? reject(error) : resolve(distanceMatrix);
    };
  });

}

export async function getCombinedDistanceMatrix(columns: DG.Column[], method: 'euclidean' | 'manhattan' = 'euclidean') {
    const promises: Promise<unknown>[] = [];
    
    columns.forEach(col => {
        if(col.type === DG.TYPE.FLOAT || col.type === DG.TYPE.INT) {
            promises.push(getDistanceMatrixForNumerics(col.toList() as unknown as Float32Array));
        }
        else //if(col.semType ==='Macromolecule') {
            {promises.push(getDistanceMatrixForSequences(col.toList() as unknown as string[]));
        }
    });

    // return the array of resolved promises values after all are resolved
    return Promise.all(promises).then((values) => {
      console.log(values);
        // if (values.length === 1)
        //     return values[0] as DistanceMatrix;
            
        // else
        //     return DistanceMatrix.combine(values as DistanceMatrix[]);
    }).catch((error) => {
        grok.shell.error(error);
    });
}