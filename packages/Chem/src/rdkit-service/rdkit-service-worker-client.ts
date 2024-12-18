//@ts-ignore
import RdKitWorkerClass from '../rdkit.worker.ts'; // .ts!
import {WORKER_CALL} from './rdkit-service-worker-api';
import {WorkerMessageBusClient} from '../worker-message-bus-client';
import {Fingerprint} from '../utils/chem-common';
import {RuleId} from '../panels/structural-alerts.js';
import {SubstructureSearchType} from '../constants.js';

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

  searchSubstructure =
    async (query: string, queryMolBlockFailover: string, molecules: string[], searchType: SubstructureSearchType) =>
      this.call(WORKER_CALL.SEARCH_SUBSTRUCTURE, [query, queryMolBlockFailover, molecules, searchType]);

  freeMoleculesStructures = async () =>
    this.call(WORKER_CALL.FREE_MOLECULES_STRUCTURES);

  getFingerprints = async (fingerprintType: Fingerprint, molecules?: string[], getCanonicalSmiles?: boolean) =>
    this.call(WORKER_CALL.GET_FINGERPRINTS, [fingerprintType, molecules, getCanonicalSmiles]);

  convertMolNotation = async (molecules: string[], targetNotation: string) =>
    this.call(WORKER_CALL.CONVERT_MOL_NOTATION, [molecules, targetNotation]);

  getStructuralAlerts = async (alerts: {[rule in RuleId]?: string[]}, molecules?: string[]) =>
    this.call(WORKER_CALL.GET_STRUCTURAL_ALERTS, [alerts, molecules]);

  rGroupAnalysis = async (molecules: string[], coreMolecule: string, coreIsQMol: boolean, options?: string) =>
    this.call(WORKER_CALL.R_GROUP_ANALYSIS, [molecules, coreMolecule, coreIsQMol, options]);

  mmpGetFragments = async (molecules: string[]) =>
    this.call(WORKER_CALL.MMP_GET_FRAGMENTS, [molecules]);

  mmpLinkFragments = async (cores: string [], fragments: string []) =>
    this.call(WORKER_CALL.MMP_LINK_FRAGMENTS, [cores, fragments]);

  mmpGetMcs = async (molecules: [string, string][]) =>
    this.call(WORKER_CALL.MMP_GET_MCS, [molecules]);

  invalidateCache = async () => this.call(WORKER_CALL.INVALIDATE_CACHE);

  setTerminateFlag = async (flag: boolean) => this.call(WORKER_CALL.SET_TERMINATE_FLAG, [flag]);

  mostCommonStructure = async (
    molecules: string[], exactAtomSearch: boolean, exactBondSearch: boolean,
  ): Promise<string> =>
    this.call(WORKER_CALL.MOST_COMMON_STRUCTURE, [molecules, exactAtomSearch, exactBondSearch]);

  getInverseSubstructuresAndAlign = async (cores: string[], from: string[], to: string[]) =>
    this.call(WORKER_CALL.INVERSE_SUBSTRUCTURE_AND_ALIGN, [cores, from, to]);
}
