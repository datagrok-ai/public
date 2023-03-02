/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// import {_package} from '../package';

const enum QUERY_NAME {
  START = 'Sequence',
  USERS = 'Users',
  ICDS = 'ICDs',
  SOURCES = 'Sources',
  SALTS = 'Salts'

}
export abstract class DataLoaderBase {
  abstract getUsers(): Promise<DG.DataFrame>;
  // abstract get users(): DG.DataFrame;
  // abstract get users(): DG.DataFrame;
  // abstract get users(): DG.DataFrame;
  // abstract get users(): DG.DataFrame;
}

/** Load data for users, salts, ICDs, IDPs and sources. In order for that to
 * work, we suppose that the queries have standard names of the form 'Sequence
 * Translator <name>' */
export class DataLoaderDB extends DataLoaderBase {
  public async init(): Promise<void> {
    const dataQueryList = (await grok.dapi.queries.list())
      .filter((dq) => dq.name.startsWith(QUERY_NAME.START));
    const usersQuery = dataQueryList.find(
      (dq) => dq.name.endsWith(QUERY_NAME.USERS)
    )!;
    this.usersDf = await this.getDatafame(usersQuery);
    if (this.usersDf === undefined)
      throw new Error('users undefined!');
  };

  private usersDf?: DG.DataFrame;

  private async getDatafame(query: DG.DataQuery) {
    const df: DG.DataFrame = (await query.prepare().call())
      .getOutputParamValue();
    console.log('inside getDataframe: df:', df.toString());
    return df;
  }

  public async getUsers(): Promise<DG.DataFrame> {
    if (this.usersDf === undefined)
      await this.init();

    return this.usersDf!;
  }
}
