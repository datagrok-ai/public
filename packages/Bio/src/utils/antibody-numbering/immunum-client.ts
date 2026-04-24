import type {
  ImmunumNumberingRow, ImmunumWorkerRequest, ImmunumWorkerResponse,
} from './types';

/** Short-lived worker: we spin it up per call and tear it down immediately
 *  after. Numbering is a one-shot batch operation — keeping the worker alive
 *  would pin the immunum WASM instance (≈700 KB) in memory indefinitely. */
function spawnWorker(): Worker {
  return new Worker(new URL('./immunum.worker', import.meta.url));
}

function callOnce(worker: Worker, req: ImmunumWorkerRequest): Promise<ImmunumWorkerResponse> {
  return new Promise((resolve, reject) => {
    const ch = new MessageChannel();
    ch.port1.onmessage = ({data}) => {
      ch.port1.close();
      resolve(data as ImmunumWorkerResponse);
    };
    ch.port1.onmessageerror = (err) => {
      ch.port1.close();
      reject(err);
    };
    worker.postMessage({req}, [ch.port2]);
  });
}

/** Runs immunum numbering on a batch of sequences inside a web worker. Spawns a
 *  fresh worker for this call and terminates it before returning so the WASM
 *  instance is freed. Throws on WASM/init errors; individual per-row errors
 *  are attached to each row's `error` field. */
export async function numberSequencesWithImmunum(
  sequences: string[],
  scheme: string,
  chains?: string[],
  minConfidence?: number | null,
): Promise<ImmunumNumberingRow[]> {
  const worker = spawnWorker();
  try {
    const resp = await callOnce(worker, {op: 'number', sequences, scheme, chains, minConfidence});
    if (!resp.ok) throw new Error(resp.error);
    return resp.rows ?? [];
  } finally {
    worker.terminate();
  }
}
