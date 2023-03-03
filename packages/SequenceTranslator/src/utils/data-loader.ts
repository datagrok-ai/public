/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// import {_package} from '../package';

const enum QUERY_NAME {
  START = 'sequence',
  USERS = 'users',
  ICDS = 'icds',
  IDPS = 'idps',
  SOURCES = 'sources',
  SALTS = 'salts'
}

export abstract class DataLoaderBase {
  abstract getUsers(): Promise<DG.DataFrame>;
  abstract getICDs(): Promise<DG.DataFrame>;
  abstract getIDPs(): Promise<DG.DataFrame>;
  abstract getSources(): Promise<DG.DataFrame>;
  abstract getSalts(): Promise<DG.DataFrame>;
}

/** Load data for users, salts, ICDs, IDPs and sources. In order for that to
 * work, we suppose that the queries have standard names of the form 'Sequence
 * Translator <name>' */
export class DataLoaderDB extends DataLoaderBase {
  public async init(): Promise<void> {
    let dataQueryList = (await grok.dapi.queries.include('params').list());
    console.log('non-fitered list:', dataQueryList.map((dq) => dq.name));
    dataQueryList = dataQueryList.filter((dq) => dq.name.toLowerCase().startsWith(QUERY_NAME.START));
    console.log('filtered list:', dataQueryList.map((dq) => dq.name));
    console.log('query list length:', dataQueryList.length);
    for (const tableName of this.tableNames) {
      console.log('tableName:', tableName);
      const query = dataQueryList.find(
        (dq) => dq.name.toLowerCase().endsWith(tableName)
      )!;
      const df = await this.getDatafameFromQuery(query);
      console.log('df:', df);
      this.dataframes.set(tableName, df);
      console.log(tableName, df.toCsv().toString());
    }
  };

  // todo: define a map of names to DFs
  private tableNames = [QUERY_NAME.USERS, QUERY_NAME.ICDS, QUERY_NAME.IDPS, QUERY_NAME.SOURCES, QUERY_NAME.SALTS];

  private dataframes = new Map<QUERY_NAME, DG.DataFrame>();

  private async getDatafameFromQuery(query: DG.DataQuery) {
    const df: DG.DataFrame = (await query.prepare().call())
      .getOutputParamValue();
    return df;
  }

  private async getDataframeByName(name: QUERY_NAME): Promise<DG.DataFrame> {
    if (this.dataframes === undefined)
      await this.init();

    return this.dataframes.get(name)!;
  }

  public async getUsers(): Promise<DG.DataFrame> {
    return this.getDataframeByName(QUERY_NAME.USERS);
  }

  public async getICDs(): Promise<DG.DataFrame> {
    return this.getDataframeByName(QUERY_NAME.ICDS);
  }

  public async getIDPs(): Promise<DG.DataFrame> {
    return this.getDataframeByName(QUERY_NAME.IDPS);
  }

  public async getSalts(): Promise<DG.DataFrame> {
    return this.getDataframeByName(QUERY_NAME.SALTS);
  }

  public async getSources(): Promise<DG.DataFrame> {
    return this.getDataframeByName(QUERY_NAME.SOURCES);
  }
}
