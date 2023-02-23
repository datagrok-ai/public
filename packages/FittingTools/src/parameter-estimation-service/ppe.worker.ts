import {WORKER_CALL} from './pe-service-worker-api';
import {PeServiceWorker as ServiceWorkerClass} from './pe-service-worker';

const ctx: Worker = self as any;
let _peServiceWorker: ServiceWorkerClass | null = null;

ctx.addEventListener('message', async (e: any) => {
  const {op, args} = e.data;
  const port = e.ports[0];
  if (op === 'module::init') {
    console.log('PE (worker) initialized');
    _peServiceWorker = new ServiceWorkerClass();
    port.postMessage({op: op, retval: null});
  } else if (op === WORKER_CALL.FIT) {
    const result = await _peServiceWorker!.fit(args[0], args[1], args[2]);
    port.postMessage({op: op, retval: result});
  }
});

