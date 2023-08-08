import RdKitWorkerClass from '../rdkit.worker.ts'; // .ts!
import {WORKER_CALL} from './rdkit-service-worker-api';
import {WorkerMessageBusClient} from '../worker-message-bus-client';
import {Fingerprint} from '../utils/chem-common';
import { RuleId } from '../panels/structural-alerts.js';

export class RdKitServiceWorkerClient extends WorkerMessageBusClient {
  constructor() {
    super(new RdKitWorkerClass());
  }

  moduleInit = async (pathToRdkit: string) =>
    this.call('module::init', [pathToRdkit]);

  /** Creates RDMols for the specified {@link molecules}.
   * They will be used for subsequent substructure search, or calculation of fingerprints.
   * Returns a number of malformed molecules. */
  initMoleculesStructures = async (molecules: string[]): Promise<number> =>
    this.call(WORKER_CALL.INIT_MOLECULES_STRUCTURES, [molecules]);

  searchSubstructure = async (query: string, queryMolBlockFailover: string, molecules: string[]) =>
    this.call(WORKER_CALL.SEARCH_SUBSTRUCTURE, [query, queryMolBlockFailover, molecules]);
    
  freeMoleculesStructures = async () =>
    this.call(WORKER_CALL.FREE_MOLECULES_STRUCTURES);

  getFingerprints = async (fingerprintType: Fingerprint, molecules?: string[], getCanonicalSmiles?: boolean) =>
    this.call(WORKER_CALL.GET_FINGERPRINTS, [fingerprintType, molecules, getCanonicalSmiles]);

  convertMolNotation = async (targetNotation: string, bitset?: boolean[]) =>
    this.call(WORKER_CALL.CONVERT_MOL_NOTATION, [targetNotation, bitset]);

  getStructuralAlerts = async (alerts: {[rule in RuleId]?: string[]}, molecules?: string[]) =>
    this.call(WORKER_CALL.GET_STRUCTURAL_ALERTS, [alerts, molecules]);

  invalidateCache = async () => this.call(WORKER_CALL.INVALIDATE_CACHE);
}
