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
    let out: DistanceMatrix | null = null;
    for(const col of columns) {
      let values: Float32Array;
      if(col.type === DG.TYPE.FLOAT || col.type === DG.TYPE.INT) {
          values = await getDistanceMatrixForNumerics(col.toList() as unknown as Float32Array);
      }
      else if(col.semType ==='Macromolecule') {
          values = await getDistanceMatrixForSequences(col.toList() as unknown as string[]);
      }
      else throw new TypeError('Unsupported column type');

      if(!out){
        out = new DistanceMatrix(values);
        if(columns.length > 1){
          out.normalize();
          out.square();
        }
      }
      else {
        let newMat: DistanceMatrix | null = new DistanceMatrix(values);
        newMat.normalize();
        switch(method) {
          case 'manhattan':
            out.add(newMat);
          default:
            newMat.square();
            out.add(newMat);
        }
        // remove reference
        newMat = null;
      }
    }
    if (method === 'euclidean'){
      out?.sqrt();
    }
    return out;
}
