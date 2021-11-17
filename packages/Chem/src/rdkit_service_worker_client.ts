import Worker from "./rdkit.worker.ts"; // .ts!
import {WorkerMessageBusClient} from './worker-message-bus-client';

export class RdKitServiceWorkerClient extends WorkerMessageBusClient {
  constructor() { super(new Worker()); }
  moduleInit = async (pathToRdkit: string) =>
    this.call('module::init', [pathToRdkit]);
  substructInit = async (dict: any) =>
    this.call('substructLibrary::init', [dict]);
  substructSearch = async (query: string, querySmarts: string) =>
    this.call('substructLibrary::search', [query, querySmarts]);
  substructDeinit = async () =>
    this.call('substructLibrary::deinit');
}