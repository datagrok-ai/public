import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {MAX_MCS_ROW_COUNT} from '../constants';


export async function getMCS(molecules: DG.Column<string>, exactAtomSearch: boolean, exactBondSearch: boolean): Promise<string> {
  if (molecules.length > MAX_MCS_ROW_COUNT) {
    grok.shell.warning(`Too many rows, maximum for MCS is ${MAX_MCS_ROW_COUNT}`);
    return '';
  }
  return new Promise((resolve, reject) => {
    const worker = new Worker(new URL('./most-common-subs-worker.ts', import.meta.url));
    // we cannot send molecules column directly to worker, because it is not a serializable object
    worker.postMessage({molecules: molecules.toList(), exactAtomSearch, exactBondSearch});
    worker.onmessage = (event) => {
      worker.terminate();
      const {error, smarts} = event.data; 
      error ? reject(error) : resolve(smarts);
    };
  });
}

