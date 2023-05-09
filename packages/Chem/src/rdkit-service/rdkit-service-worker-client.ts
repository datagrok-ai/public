import RdKitWorkerClass from '../rdkit.worker.ts'; // .ts!
import {WORKER_CALL} from './rdkit-service-worker-api';
import {WorkerMessageBusClient} from '../worker-message-bus-client';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {Fingerprint} from '../utils/chem-common';
import {RDModule} from "@datagrok-libraries/chem-meta/src/rdkit-api";

export class RdKitServiceWorkerClient extends WorkerMessageBusClient {
  constructor() {
    super(new RdKitWorkerClass());
  }

  moduleInit = async (pathToRdkit: string) =>
    this.call('module::init', [pathToRdkit]);

  /** Creates RDMols for the specified {@link molecules}.
   * They will be used for subsequent substructure search, or calculation of fingerprints.
   * Returns a number of malformed molecules. */
  initMoleculesStructures = async (dict: string[]): Promise<number> =>
    this.call(WORKER_CALL.INIT_MOLECULES_STRUCTURES, [dict]);

  searchSubstructure = async (query: string, queryMolBlockFailover: string, bitset?: boolean[]) =>
    this.call(WORKER_CALL.SEARCH_SUBSTRUCTURE, [query, queryMolBlockFailover, bitset]);

  freeMoleculesStructures = async () =>
    this.call(WORKER_CALL.FREE_MOLECULES_STRUCTURES);

  getFingerprints = async (fingerprintType: Fingerprint) =>
    this.call(WORKER_CALL.GET_FINGERPRINTS, [fingerprintType]);

  convertMolNotation = async (targetNotation: string, bitset?: boolean[]) =>
    this.call(WORKER_CALL.CONVERT_MOL_NOTATION, [targetNotation, bitset]);
}
