/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import {_package} from '../../package';
import {Fingerprint, defaultMorganFpLength} from '../../utils/chem-common';
import {getRdKitService} from '../../utils/chem-common-rdkit';

/** Run BitBIRCH clustering on a molecule column. Returns an int column with cluster assignments. */
export async function bitbirchWorker(
  molCol: DG.Column,
  threshold: number = 0.65,
  fingerprintType: Fingerprint = Fingerprint.Morgan,
  nFeatures: number = defaultMorganFpLength,
  branchingFactor: number = 50,
): Promise<DG.Column> {
  const nBytesPerFp = Math.ceil(nFeatures / 8);
  const rowCount = molCol.length;
  const rdkitService = await getRdKitService();
  const fps = (await rdkitService.getFingerprints(fingerprintType, molCol.toList(), false)).fps;
  // Compute fingerprints and pack into contiguous buffer
  const packedFps = new Uint8Array(rowCount * nBytesPerFp);
  let validCount = 0;
  // Map from position in valid-only array back to original row index
  const validIndices: number[] = [];

  for (let i = 0; i < fps.length; i++) {
    try {
      if (!fps[i]) continue;
      packedFps.set(fps[i]!, validCount * nBytesPerFp);
      validIndices.push(i);
      validCount++;
    } catch (_) {
      continue;
    }
  }

  if (validCount === 0) {
    const col = DG.Column.int('Cluster (BitBIRCH)', rowCount);
    col.init((_i) => DG.INT_NULL);
    return col;
  }

  // Run BitBIRCH on valid fingerprints only
  const validFps = packedFps.slice(0, validCount * nBytesPerFp);
  const worker = new Worker(new URL('./bitbirch-worker.ts', import.meta.url));
  return new Promise((resolve, reject) => {
    worker.onmessage = (event: MessageEvent) => {
      const {assignments, error}: {assignments?: Int32Array, error?: string} = event.data;
      try {
        if (error || !assignments) {
          worker.terminate();
          reject(new Error(`BitBIRCH worker error: ${error}`));
          return;
        }
        const result = DG.Column.int('Cluster (BitBIRCH)', rowCount);
        result.init((_i) => DG.INT_NULL);
        for (let j = 0; j < validCount; j++)
          result.set(validIndices[j], assignments[j]);
        worker.terminate();
        resolve(result);
      } catch (e) {
        worker.terminate();
        reject(e instanceof Error ? e : new Error(String(e)));
      }
    };
    worker.onerror = (err) => {
      worker.terminate();
      reject(new Error(`BitBIRCH worker error: ${err.message}`));
    };
    worker.postMessage({threshold, branchingFactor, nFeatures, fps: validFps, fpCount: validCount, webRoot: _package.webRoot});
  });
}
