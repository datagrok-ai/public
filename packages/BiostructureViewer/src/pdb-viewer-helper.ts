type Tracks = {[kind: string]: number[]};
type _Tracks = {[kind: string]: Set<number>};

/**
 * Polymer entity instance (a.k.a chain) data.
 *
 * @interface Instance
 */
interface Instance {
  id: string,
  tracks: Tracks;
}

/**
 * Polymer entity data.
 *
 * @interface Entity
 */
interface Entity {
  id: string;
  sequence: string;
  instances: Instance[];
}

/**
 * RCSB REST API URLs.
 *
 * @enum {number}
 */
enum RCSBRESTAPI {
  file = 'https://files.rcsb.org/download/{entry_id}.pdb',
  entry = 'https://data.rcsb.org/rest/v1/core/entry/{entry_id}',
  entity = 'https://data.rcsb.org/rest/v1/core/polymer_entity/{entry_id}/{entity_id}',
  instance = 'https://data.rcsb.org/rest/v1/core/polymer_entity_instance/{entry_id}/{asym_id}',
}

/**
 * Can be used to retrieve primary, secondary, and tertiary structure
 * of polymers from RCSB PDB. Primarily designed for proteins.
 *
 * @export
 * @class PDBViewerHelper
 */
export class PDBViewerHelper {
  protected RCSB = RCSBRESTAPI;
  protected secondaryKinds = ['SHEET', 'HELIX_P'];
  protected pdbID: string;
  protected pdbBody: string;
  protected items: Entity[];

  /**
   * Creates an instance of PDBViewerHelper.
   * @param {string} pdbID PDB ID of a polymer.
   * @memberof PDBViewerHelper
   */
  constructor(pdbID: string) {
    this.pdbID = pdbID.toLocaleLowerCase();
    this.pdbBody = '';
    this.items = [];
  }

  /**
   * Fetches all pieces of data.
   *
   * @memberof PDBViewerHelper
   */
  async fetchInfo() {
    this.pdbBody = await this._getPDBAsText();

    const entry = await this._getPDBEntry();
    const entityIDs = entry['rcsb_entry_container_identifiers']['polymer_entity_ids'];

    this.items = await this._cycleThroughPDBEntities(entityIDs);
  }

  /**
   * Downloads structure file as a string.
   *
   * @protected
   * @return {Promise<string>} String containing structure file.
   * @memberof PDBViewerHelper
   */
  protected async _getPDBAsText(): Promise<string> {
    const link = this.RCSB.file.replace('{entry_id}', this.pdbID);
    return await _fetchTextFromURL(link);
  }

  /**
   * Downloads information on PDB entry.
   *
   * @protected
   * @return {Promise<any>} PDB entry object.
   * @memberof PDBViewerHelper
   */
  protected async _getPDBEntry(): Promise<any> {
    const link = this.RCSB.entry.replace('{entry_id}', this.pdbID);
    return await _fetchTextFromURL(link, true);
  }

  /**
   * Downloads PDB entity information.
   *
   * @protected
   * @param {string} entityID ID of an entity.
   * @return {Promise<any>} Entity object.
   * @memberof PDBViewerHelper
   */
  protected async _getPDBEntity(entityID: string): Promise<any> {
    const link = this.RCSB.entity
      .replace('{entry_id}', this.pdbID)
      .replace('{entity_id}', entityID);
    return await _fetchTextFromURL(link, true);
  }

  /**
   * Iterates over entities to obtain primary and secondary structure.
   *
   * @protected
   * @param {string[]} entityIDs Entities' ID.
   * @return {Promise<Entity[]>} Entities with polymer sequence and secondary structure.
   * @memberof PDBViewerHelper
   */
  protected async _cycleThroughPDBEntities(entityIDs: string[]): Promise<Entity[]> {
    const entities: Entity[] = [];

    for (const entityID of entityIDs) {
      const entity = await this._getPDBEntity(entityID);

      const asymIDs = entity['rcsb_polymer_entity_container_identifiers']['asym_ids'];
      const seq = entity['entity_poly']['pdbx_seq_one_letter_code'];
      const instances: Instance[] = [];

      for (const asymID of asymIDs) {
        const tracks = await this._filterPDBPolymerEntityInstances(asymID);
        instances.push({id: asymID, tracks: tracks});
      }

      entities.push({id: entityID, sequence: seq, instances: instances});
    }
    return entities;
  }

  /**
   * Filters polymer chains instance to get sheet and helix tracks.
   *
   * @protected
   * @param {string} asymID Chain ID of an entity.
   * @return {Promise<Tracks>} Positions with secondary structure information assigned.
   * @memberof PDBViewerHelper
   */
  protected async _filterPDBPolymerEntityInstances(asymID: string): Promise<Tracks> {
    const link = this.RCSB.instance
      .replace('{entry_id}', this.pdbID)
      .replace('{asym_id}', asymID);
    const instance = await _fetchTextFromURL(link, true);

    const tracks: Tracks = {};
    const _tracks: _Tracks = {};

    for (const annot of instance['rcsb_polymer_instance_feature']) {
      const annotType = annot['type'];

      if (this.secondaryKinds.includes(annotType)) {
        _tracks[annotType] = _tracks[annotType] ? _tracks[annotType] : new Set();

        for (const pos of annot['feature_positions']) {
          const range = _genRange(pos['beg_seq_id'], pos['end_seq_id']);

          for (const i of range) {
            _tracks[annotType].add(i);
          }
        }
      }
    }

    for (const [k, v] of Object.entries(_tracks)) {
      tracks[k] = Array.from(v.values()).sort((a: number, b: number) => a - b);
    }
    return tracks;
  }

  /**
   * Tertiary structure as a string.
   *
   * @readonly
   * @type {string}
   * @memberof PDBViewerHelper
   */
  get body(): string {
    if (this.pdbBody.length == 0) {
      throw new Error('Call fetchInfo before getting file body.');
    }
    return this.pdbBody;
  }

  /**
   * Primary and secondary structure.
   *
   * @readonly
   * @type {Entity[]}
   * @memberof PDBViewerHelper
   */
  get entities(): Entity[] {
    return this.items;
  }
}

/**
 * Helper function to request information by URL.
 *
 * @param {string} url
 * @param {boolean} [asJSON=false]
 * @return {*}  {(Promise<string | any>)}
 */
async function _fetchTextFromURL(url: string, asJSON = false): Promise<string | any> {
  const response = await fetch(url);
  return await (asJSON ? response.json() : response.text());
}

/**
 * Generates array from a range.
 *
 * @param {number} begin First position to include.
 * @param {number} end Last position to include.
 * @return {number[]} Array corresponding the range given.
 */
function _genRange(begin: number, end: number): number[] {
  const nItems = end - begin;
  return new Array(nItems).fill(0).map((_, i) => begin + i);
}
