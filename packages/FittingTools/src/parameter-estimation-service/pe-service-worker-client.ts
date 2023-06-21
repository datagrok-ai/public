import {WORKER_CALL} from './pe-service-worker-api';
import {WorkerMessageBusClient} from './worker-message-bus-client';
import PeWorkerClass from './ppe.worker.ts'; // .ts!

export class PeServiceWorkerClient extends WorkerMessageBusClient {
  constructor() {
    super(new PeWorkerClass());
  }

  moduleInit = async (pathToRdkit: string) =>
    this.call('module::init', [pathToRdkit]);

  fit = async (x: number[], y: number[], counts: number[]) =>
    this.call(WORKER_CALL.FIT, [x, y, counts]);
}
