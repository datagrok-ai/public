import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import Dexie from 'dexie';
import {keys} from 'ts-transformer-keys';


export interface ICached {

}

export interface IListVersion {
  list_id: number;
  name: string;
  version: string;
}

export interface IDbVersion extends IListVersion {

}

export interface IObject extends ICached {
  id?: number;
  key: string;
  value: any;
}

export class MlbDatabase extends Dexie {
  private readonly serverListVersionDf!: DG.DataFrame;

  get enabled(): boolean { return !!this.serverListVersionDf;}

  list_version!: Dexie.Table<IListVersion, number>;
  object!: Dexie.Table<IObject, number>;

  constructor(serverListVersionDf: DG.DataFrame) {
    super('MlbDatabase');
    this.serverListVersionDf = serverListVersionDf;
    if (this.enabled) {
      this.version(1).stores({
        list_version: 'list_id, name, version',

        scheme: 'scheme_id, scheme',
        cdr: 'cdr_id, cdr',
        antigen: 'id, antigen, antigen_ncbi_id, antigen_gene_symbol',
        vid: 'v_id',
        vidObsPtm: 'v_id',

        object: '++id, key, value',
      });

      this.version(2).stores({
        list_version: 'list_id, name, version',
        object: '++id, key, value',

        scheme: null, cdr: null, antigen: null, vid: null, vidObsPtm: null,
      }).upgrade(async () => {
        for (const key in ['scheme', 'cdr', 'antigen', 'vid', 'v_id']) {
          const cacheVersion = (await this.list_version.where({'name': key}).first());
          await this.list_version.delete(cacheVersion.list_id);
        }
      });

    }
  }

  async getObject<TData>(key: string, getServerObject: () => Promise<TData>): Promise<TData> {
    const dbVersion: IDbVersion = this.enabled ? this.getServerVersion(key) : null;

    // read version of cached list if cache is enabled
    const cacheVersion = this.enabled ? (await this.list_version.where({'name': key}).first())?.version : null;

    if (!this.enabled || cacheVersion != dbVersion.version) {
      const serverObj: TData = await getServerObject();

      if (this.enabled) {
        await Promise.all([
          this.object.put({key: key, value: serverObj}),
          this.list_version.put(dbVersion, dbVersion.list_id)
        ]);
      }

      return serverObj;
    } else {
      const cacheObj: IObject = await this.object.where({key: key}).first();
      return cacheObj.value;
    }
  }

  async getDataFrame(key: string, getServerObject: () => Promise<DG.DataFrame>): Promise<DG.DataFrame> {
    const dbVersion: IDbVersion = this.enabled ? this.getServerVersion(key) : null;

    // read version of cached list if cache is enabled
    const cacheVersion = this.enabled ? (await this.list_version.where({'name': key}).first())?.version : null;

    if (!this.enabled || cacheVersion != dbVersion.version) {
      const serverDf: DG.DataFrame = await getServerObject();

      if (this.enabled) {
        const serverDfBa: Uint8Array = serverDf.toByteArray();
        await Promise.all([
          this.object.put({key: key, value: serverDfBa}),
          this.list_version.put(dbVersion, dbVersion.list_id)
        ]);
      }

      return serverDf;
    } else {
      const cacheObj: IObject = await this.object.where({key: key}).first();
      const cacheDfBa: Uint8Array = cacheObj.value;
      const cacheDf: DG.DataFrame = DG.DataFrame.fromByteArray(cacheDfBa);
      return cacheDf;
    }
  }

  getServerVersion(name: string): IDbVersion {
    const df: DG.DataFrame = this.serverListVersionDf.rows.match({name: name}).toDataFrame();
    if (df.rowCount == 0)
      throw new Error(`Name '${name}' not found in server side versions.`);
    return Object.assign({},
      ...(['list_id', 'name', 'version'].map((fn) => ({[fn]: df.get(fn, 0)}))));
  }
}
