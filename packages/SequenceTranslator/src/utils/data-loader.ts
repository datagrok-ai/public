/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

const enum QUERY_NAME {
  START = 'sequence',
  USERS = 'users',
  ICDS = 'icds',
  IDPS = 'idps',
  SOURCES = 'sources',
  SALTS = 'salts'
}

export abstract class DataLoaderBase {
  abstract get Users(): DG.DataFrame;

  abstract get ICDs(): DG.DataFrame;

  abstract get IDPs(): DG.DataFrame;

  abstract get Sources(): DG.DataFrame;

  abstract get Salts(): DG.DataFrame;
}

/** Load data for entities as Users, Salts, ICDs, IDPs and Sources.
 * In order for that to work, we suppose that the queries have tags
 * 'app-SequenceTranslator' and 'entity-<Entity>' (e.g. 'entity-Users').
 */
export class DataLoaderDB extends DataLoaderBase {
  private _users: DG.DataFrame;
  private _ICDs: DG.DataFrame;
  private _IDPs: DG.DataFrame;
  private _sources: DG.DataFrame;
  private _salts: DG.DataFrame;

  get Users(): DG.DataFrame { return this._users; }

  get ICDs(): DG.DataFrame { return this._ICDs; }

  get IDPs(): DG.DataFrame { return this._IDPs; }

  get Sources(): DG.DataFrame { return this._sources; }

  get Salts(): DG.DataFrame { return this._salts; }

  public async init(): Promise<void> {
    // Server side filter by tag requires '#' sign prefix
    const dataQueryList: DG.DataQuery[] = (await grok.dapi.queries.include('connection,params,entityTags')
      .filter('#app-SequenceTranslator').list());

    const usersDQ: DG.DataQuery = dataQueryList.find((dq) => dq.hasTag('entity-Users'))!;
    const icdsDQ: DG.DataQuery = dataQueryList.find((dq) => dq.hasTag('entity-ICDs'))!;
    const idpsDQ: DG.DataQuery = dataQueryList.find((dq) => dq.hasTag('entity-IDPs'))!;
    const sourcesDQ: DG.DataQuery = dataQueryList.find((dq) => dq.hasTag('entity-Sources'))!;
    const saltsDQ: DG.DataQuery = dataQueryList.find((dq) => dq.hasTag('entity-Salts'))!;

    await Promise.all([
      usersDQ.prepare().call().then((fc) => { this._users = fc.getOutputParamValue(); }),
      icdsDQ.prepare().call().then((fc) => { this._ICDs = fc.getOutputParamValue(); }),
      idpsDQ.prepare().call().then((fc) => { this._IDPs = fc.getOutputParamValue(); }),
      sourcesDQ.prepare().call().then((fc) => {this._sources = fc.getOutputParamValue(); }),
      saltsDQ.prepare().call().then((fc) => { this._salts = fc.getOutputParamValue(); }),
    ]);
  };
}
