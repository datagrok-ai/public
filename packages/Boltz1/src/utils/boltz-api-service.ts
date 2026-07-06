import * as grok from 'datagrok-api/grok';

import Boltz from 'boltz-api';
import type {Fetch} from 'boltz-api/internal/builtin-types';
import type {StructureAndBindingStartParams} from 'boltz-api/resources/predictions/structure-and-binding';
import type {AdmeStartParams} from 'boltz-api/resources/predictions/adme';
import type {DesignStartParams as SmallMoleculeDesignStartParams} from 'boltz-api/resources/small-molecule/design';
import type {
  LibraryScreenStartParams as SmallMoleculeLibraryScreenStartParams,
} from 'boltz-api/resources/small-molecule/library-screen';
import type {DesignStartParams as ProteinDesignStartParams} from 'boltz-api/resources/protein/design';
import type {
  LibraryScreenStartParams as ProteinLibraryScreenStartParams,
} from 'boltz-api/resources/protein/library-screen';

import {_package} from '../package';
import {
  BOLTZ_API_BASE_URL, BOLTZ_API_KEY_PARAM, BOLTZ_POLL_INTERVAL_MS, BoltzJob, TERMINAL_STATUSES,
} from './boltz-api-constants';

const proxyFetch: Fetch = (input, init) => {
  const url = input instanceof Request ? input.url : input.toString();
  const headers: Record<string, string> = {};
  new Headers(init?.headers).forEach((value, key) => headers[key] = value);
  return grok.dapi.fetchProxy(url, {...init, headers});
};

export async function getInitializedService(): Promise<BoltzService> {
  const credentials = await _package.getCredentials();
  const apiKey = credentials?.parameters[BOLTZ_API_KEY_PARAM];
  const svc = BoltzService.getInstance();
  await svc.init(apiKey);
  return svc;
}

/** Singleton over the Boltz hosted API: caches the authenticated client and runs each task end-to-end. */
export class BoltzService {
  private static _instance: BoltzService | null = null;
  private boltzClient: Boltz | undefined;
  private apiKey: string | undefined;

  private constructor() {}

  static getInstance(): BoltzService {
    return BoltzService._instance ??= new BoltzService();
  }

  async init(apiKey: string | undefined): Promise<void> {
    if (this.boltzClient)
      return;
    if (!apiKey)
      throw new Error('Boltz API key is not set in the Boltz1 package credentials.');
    this.apiKey = apiKey;
    this.boltzClient = new Boltz({apiKey, baseURL: BOLTZ_API_BASE_URL, fetch: proxyFetch});
  }

  private get client(): Boltz {
    if (!this.boltzClient)
      throw new Error('BoltzService is not initialized; call init() first.');
    return this.boltzClient;
  }

  /** Single predictions (structure & binding, ADME): start, then poll until the inline output is ready. */
  async predictStructureAndBinding(request: StructureAndBindingStartParams) {
    const {id} = await this.client.predictions.structureAndBinding.start(request);
    return this.pollJob(() => this.client.predictions.structureAndBinding.retrieve(id));
  }

  async predictAdme(request: AdmeStartParams) {
    const {id} = await this.client.predictions.adme.start(request);
    return this.pollJob(() => this.client.predictions.adme.retrieve(id));
  }

  /** Pipeline runs (design & library screen): start, poll until terminal, then collect all scored results. */
  async designSmallMolecules(request: SmallMoleculeDesignStartParams) {
    const {id} = await this.client.smallMolecule.design.start(request);
    await this.pollJob(() => this.client.smallMolecule.design.retrieve(id));
    return this.collectResults(this.client.smallMolecule.design.listResults(id));
  }

  async screenSmallMoleculeLibrary(request: SmallMoleculeLibraryScreenStartParams) {
    const {id} = await this.client.smallMolecule.libraryScreen.start(request);
    await this.pollJob(() => this.client.smallMolecule.libraryScreen.retrieve(id));
    return this.collectResults(this.client.smallMolecule.libraryScreen.listResults(id));
  }

  async designProteins(request: ProteinDesignStartParams) {
    const {id} = await this.client.protein.design.start(request);
    await this.pollJob(() => this.client.protein.design.retrieve(id));
    return this.collectResults(this.client.protein.design.listResults(id));
  }

  async screenProteinLibrary(request: ProteinLibraryScreenStartParams) {
    const {id} = await this.client.protein.libraryScreen.start(request);
    await this.pollJob(() => this.client.protein.libraryScreen.retrieve(id));
    return this.collectResults(this.client.protein.libraryScreen.listResults(id));
  }

  /** Polls a job via its retrieve method until it reaches a terminal status; throws if it failed. */
  private async pollJob<T extends BoltzJob>(
    retrieve: () => Promise<T>,
    intervalMs = BOLTZ_POLL_INTERVAL_MS,
  ): Promise<T> {
    let job = await retrieve();
    while (!TERMINAL_STATUSES.has(job.status)) {
      await new Promise((resolve) => setTimeout(resolve, intervalMs));
      job = await retrieve();
    }
    if (job.status === 'failed')
      throw new Error(`Boltz job ${job.id} failed: ${job.error?.message ?? 'unknown error'}`);
    return job;
  }

  /** Collects every page of a cursor-paginated result set into a single array. */
  private async collectResults<T>(pages: AsyncIterable<T>): Promise<T[]> {
    const results: T[] = [];
    for await (const item of pages)
      results.push(item);
    return results;
  }
}
