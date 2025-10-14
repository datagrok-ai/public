import {genRange} from '@datagrok-libraries/utils/src/vector-operations';
import {RcsbGraphQLAdapter} from './utils/rcsb-gql-adapter';
import * as DG from 'datagrok-api/dg';

type Tracks = {[kind: string]: number[]};
type _Tracks = {[kind: string]: Set<number>};

/** Polymer entity instance (a.k.a chain) data. */
interface Chain {
  id: string,
  tracks: Tracks;
}

/** Polymer entity data. */
interface Entity {
  id: string;
  sequence: string;
  chains: Chain[];
}

/** RCSB REST API URLs. */
enum RCSBRESTAPI {
  file = 'https://files.rcsb.org/download/{entry_id}.pdb',
  entry = 'https://data.rcsb.org/rest/v1/core/entry/{entry_id}',
  entity = 'https://data.rcsb.org/rest/v1/core/polymer_entity/{entry_id}/{entity_id}',
  instance = 'https://data.rcsb.org/rest/v1/core/polymer_entity_instance/{entry_id}/{asym_id}',
}

/**
 * Enhanced PdbEntry class with GraphQL capabilities.
 * Can be used to retrieve primary, secondary, and tertiary structure
 * of polymers from RCSB PDB using both REST and GraphQL APIs.
 */
export class PdbEntry {
  protected RCSB = RCSBRESTAPI;
  protected secondaryKinds = ['SHEET', 'HELIX_P'];
  protected pdbID: string;
  protected pdbBody: string;
  protected items: Entity[];
  protected useGraphQL: boolean;

  constructor(pdbID: string, useGraphQL: boolean = true) {
    this.pdbID = pdbID.toLowerCase();
    this.pdbBody = '';
    this.items = [];
    this.useGraphQL = useGraphQL;
  }

  /**
   * Fetches all pieces of data using the preferred method (GraphQL or REST).
   */
  async fetchInfo(): Promise<void> {
    if (this.useGraphQL)
      await this.fetchInfoGraphQL();
    else
      await this.fetchInfoREST();
  }

  /**
   * Fetches polymer sequences using GraphQL.
   * Much more efficient than downloading the entire structure file.
   */
  async fetchPolymerSequences(): Promise<{[entityId: string]: {sequence: string; type: string; chains: string[]}}> {
    return await RcsbGraphQLAdapter.getPolymerSequences(this.pdbID);
  }

  /**
   * Get only protein sequences using GraphQL.
   * Filters out DNA/RNA and returns a simple chain -> sequence mapping.
   */
  async getProteinSequences(): Promise<{[chainId: string]: string}> {
    return await RcsbGraphQLAdapter.getProteinSequences(this.pdbID);
  }

  /**
   * Get basic entry information using GraphQL.
   */
  async getEntryInfo(): Promise<{id: string; title?: string; resolution?: number; experimentalMethod?: string}> {
    return await RcsbGraphQLAdapter.getEntryInfo(this.pdbID);
  }

  /**
   * Create a DataFrame with sequence data using GraphQL.
   */
  async createSequenceDataFrame(): Promise<DG.DataFrame> {
    return await RcsbGraphQLAdapter.createSequenceDataFrame(this.pdbID);
  }


  /**
   * Enhanced fetchInfo using GraphQL for polymer data and REST for structure file.
   * This hybrid approach gets sequences efficiently while still providing full structure access.
   */
  private async fetchInfoGraphQL(): Promise<void> {
    // Get basic info and sequences via GraphQL
    const [entryInfo, polymerData] = await Promise.all([
      this.getEntryInfo(),
      this.fetchPolymerSequences()
    ]);

    // Convert GraphQL data to the existing Entity format
    this.items = await this.convertGraphQLToEntities(polymerData);

    // Optionally fetch the structure file (only if needed)
    // this.pdbBody = await this._getPDBAsText();
  }

  /**
   * Original REST-based fetchInfo method (preserved for compatibility).
   */
  private async fetchInfoREST(): Promise<void> {
    this.pdbBody = await this._getPDBAsText();

    const entry = await this._getPDBEntry();
    const entityIDs = entry['rcsb_entry_container_identifiers']['polymer_entity_ids'];

    this.items = await this._cycleThroughPDBEntities(entityIDs);
  }

  /**
   * Convert GraphQL polymer data to the existing Entity interface format.
   */
  private async convertGraphQLToEntities(polymerData: {[entityId: string]: {sequence: string; type: string; chains: string[]}}): Promise<Entity[]> {
    const entities: Entity[] = [];

    for (const [entityId, data] of Object.entries(polymerData)) {
      const chains: Chain[] = [];

      // For each chain, create a Chain object
      // Note: GraphQL doesn't provide secondary structure info, so tracks will be empty
      // You could extend this to make additional GraphQL queries for secondary structure if needed
      for (const chainId of data.chains) {
        chains.push({
          id: chainId,
          tracks: {} // Empty for now - could be populated with additional GraphQL queries
        });
      }

      entities.push({
        id: entityId,
        sequence: data.sequence,
        chains: chains
      });
    }

    return entities;
  }

  /**
   * Fetch secondary structure information for a specific chain using REST.
   * This can be used to populate tracks when full secondary structure info is needed.
   */
  async fetchSecondaryStructure(asymId: string): Promise<Tracks> {
    return await this._filterPDBPolymerEntityInstances(asymId);
  }

  /** Downloads structure file as a string. */
  protected async _getPDBAsText(): Promise<string> {
    const link = this.RCSB.file.replace('{entry_id}', this.pdbID);
    return await _fetchTextFromURL(link);
  }

  /** Downloads information on the PDB entry and returns a PDB entry object. */
  protected async _getPDBEntry(): Promise<any> {
    const link = this.RCSB.entry.replace('{entry_id}', this.pdbID);
    return await _fetchTextFromURL(link, true);
  }

  /**
   * Downloads PDB entity information and returns Entity object.
   * @param {string} entityID ID of an entity.
   */
  protected async _getPDBEntity(entityID: string): Promise<any> {
    const link = this.RCSB.entity
      .replace('{entry_id}', this.pdbID)
      .replace('{entity_id}', entityID);
    return await _fetchTextFromURL(link, true);
  }

  /**
   * Iterates over entities and returns primary and secondary structure.
   * @param {string[]} entityIDs Entities' ID.
   */
  protected async _cycleThroughPDBEntities(entityIDs: string[]): Promise<Entity[]> {
    const entities: Entity[] = [];

    for (const entityID of entityIDs) {
      const entity = await this._getPDBEntity(entityID);

      const asymIDs = entity['rcsb_polymer_entity_container_identifiers']['asym_ids'];
      const seq = entity['entity_poly']['pdbx_seq_one_letter_code'];
      const chains: Chain[] = [];

      for (const asymID of asymIDs) {
        const tracks = await this._filterPDBPolymerEntityInstances(asymID);
        chains.push({id: asymID, tracks: tracks});
      }

      entities.push({id: entityID, sequence: seq, chains: chains});
    }
    return entities;
  }

  /**
   * Filters polymer chains and returns positions with secondary structure information assigned.
   * @param {string} asymID Chain ID of an entity.
   */
  protected async _filterPDBPolymerEntityInstances(asymID: string): Promise<Tracks> {
    const link = this.RCSB.instance
      .replace('{entry_id}', this.pdbID)
      .replace('{asym_id}', asymID);
    const instance = await _fetchTextFromURL(link, true);

    const tracks: Tracks = {}; // Original tracks.
    const _tracks: _Tracks = {}; // Tracks with container as a Set to collect unique numbers.

    for (const annot of instance['rcsb_polymer_instance_feature']) {
      const annotType = annot['type'];

      if (this.secondaryKinds.includes(annotType)) {
        if (_tracks[annotType] === undefined)
          _tracks[annotType] = new Set();

        for (const pos of annot['feature_positions']) {
          const range = genRange(pos['beg_seq_id'], pos['end_seq_id']);

          for (const i of range)
            _tracks[annotType].add(i);
        }
      }
    }

    // Format tracks as sorted numbers.
    for (const [k, v] of Object.entries(_tracks))
      tracks[k] = Array.from(v.values()).sort((a: number, b: number) => a - b);

    return tracks;
  }

  /** Tertiary structure as a string. */
  get body(): string {
    if (this.pdbBody.length == 0) {
      if (this.useGraphQL)
        throw new Error('Structure file not loaded. Call fetchInfoREST() or _getPDBAsText() to load structure file.');
      else
        throw new Error('Call fetchInfo() before getting file body.');
    }

    return this.pdbBody;
  }

  /** Primary and secondary structure. */
  get entities(): Entity[] {
    if (this.items.length == 0)
      throw new Error('Call fetchInfo() before getting entities.');

    return this.items;
  }

  /** Get sequences only (lightweight access). */
  get sequences(): {[entityId: string]: string} {
    const result: {[entityId: string]: string} = {};
    for (const entity of this.entities)
      result[entity.id] = entity.sequence;

    return result;
  }

  /** Switch between GraphQL and REST mode. */
  setGraphQLMode(useGraphQL: boolean): void {
    this.useGraphQL = useGraphQL;
  }

  //TODO: remove - demo use only
  set sbody(pdb: string) {
    this.pdbBody = pdb;
  }
}

/**
 * Helper function to request information by URL.
 *
 * @param {string} url URL to request.
 * @param {boolean} [asJSON=false] Whether to reply with JSON or simple string.
 * @return {(Promise<string | any>)} Result of the request.
 */
async function _fetchTextFromURL(url: string, asJSON = false): Promise<string | any> {
  const response = await fetch(url);
  return await (asJSON ? response.json() : response.text());
}
