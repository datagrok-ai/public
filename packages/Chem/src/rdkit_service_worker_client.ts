import RdKitWorkerClass from "./rdkit.worker.ts"; // .ts!
import {WorkerMessageBusClient} from './worker-message-bus-client';

export class RdKitServiceWorkerClient extends WorkerMessageBusClient {
  constructor () { super(new RdKitWorkerClass()); }
  moduleInit = async (pathToRdkit: string) =>
    this.call('module::init', [pathToRdkit]);
  substructInit = async (dict: any) =>
    this.call('substructLibrary::init', [dict]);
  substructSearch = async (query: string, querySmarts: string) =>
    this.call('substructLibrary::search', [query, querySmarts]);
  substructDeinit = async () =>
    this.call('substructLibrary::deinit');
  structuralAlertsInit = async (smarts: string[]) =>
    this.call('structuralAlerts::init', [smarts]);
  structuralAlertsGet = async (smiles: string) =>
    this.call('structuralAlerts::get', [smiles]);
}