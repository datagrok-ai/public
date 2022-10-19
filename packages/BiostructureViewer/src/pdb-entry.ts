import {genRange} from '@datagrok-libraries/utils/src/vector-operations';

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
// eslint-disable-next-line no-unused-vars
enum RCSBRESTAPI {
  // eslint-disable-next-line no-unused-vars
  file = 'https://files.rcsb.org/download/{entry_id}.pdb',
  // eslint-disable-next-line no-unused-vars
  entry = 'https://data.rcsb.org/rest/v1/core/entry/{entry_id}',
  // eslint-disable-next-line no-unused-vars
  entity = 'https://data.rcsb.org/rest/v1/core/polymer_entity/{entry_id}/{entity_id}',
  // eslint-disable-next-line no-unused-vars
  instance = 'https://data.rcsb.org/rest/v1/core/polymer_entity_instance/{entry_id}/{asym_id}',
}

/**
 * Can be used to retrieve primary, secondary, and tertiary structure
 * of polymers from RCSB PDB. Primarily designed for proteins.
 */
export class PdbEntry {
  protected RCSB = RCSBRESTAPI;
  protected secondaryKinds = ['SHEET', 'HELIX_P'];
  protected pdbID: string;
  protected pdbBody: string;
  protected items: Entity[];

  constructor(pdbID: string) {
    this.pdbID = pdbID.toLocaleLowerCase();
    this.pdbBody = '';
    this.items = [];
  }

  /** Fetches all pieces of data. */
  async fetchInfo() {
    this.pdbBody = await this._getPDBAsText();

    const entry = await this._getPDBEntry();
    const entityIDs = entry['rcsb_entry_container_identifiers']['polymer_entity_ids'];

    this.items = await this._cycleThroughPDBEntities(entityIDs);
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
   * Filters polymer chains and returns hositions with secondary structure information assigned.
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
        if (_tracks[annotType] === undefined) {
          _tracks[annotType] = new Set();
        }

        for (const pos of annot['feature_positions']) {
          const range = genRange(pos['beg_seq_id'], pos['end_seq_id']);

          for (const i of range) {
            _tracks[annotType].add(i);
          }
        }
      }
    }

    // Format tracks as sorted numbers.
    for (const [k, v] of Object.entries(_tracks)) {
      tracks[k] = Array.from(v.values()).sort((a: number, b: number) => a - b);
    }
    return tracks;
  }

  /** Tertiary structure as a string. */
  get body(): string {
    if (this.pdbBody.length == 0) {
      throw new Error('Call fetchInfo before getting file body.');
    }
    return this.pdbBody;
  }

  /** Primary and secondary structure. */
  get entities(): Entity[] {
    if (this.items.length == 0) {
      throw new Error('Call fetchInfo before getting entities.');
    }
    return this.items;
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
 * @param {boolean} [asJSON=false] Whether to reply with JSIN or simple string.
 * @return {(Promise<string | any>)} Result of the request.
 */
async function _fetchTextFromURL(url: string, asJSON = false): Promise<string | any> {
  const response = await fetch(url);
  return await (asJSON ? response.json() : response.text());
}
