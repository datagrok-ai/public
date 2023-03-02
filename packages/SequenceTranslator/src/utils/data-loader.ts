/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// import {_package} from '../package';

export abstract class DataLoaderBase {
  abstract getUsers(): Promise<DG.DataFrame>;
  // abstract get users(): DG.DataFrame;
  // abstract get users(): DG.DataFrame;
  // abstract get users(): DG.DataFrame;
  // abstract get users(): DG.DataFrame;
}

/** Load data for users, salts, ICDs, IDPs and sources. In order for that to
 * wrok, we suppose that the queries have standard names of the form "Sequence
 * Translator " */
export class DataLoaderDB extends DataLoaderBase {
  private async getDatafame(query: DG.DataQuery) {
    const c = await query.prepare().call();
    console.log('inside getDataframe: c:', c.toString());
    const df: DG.DataFrame = c.getOutputParamValue();
    console.log('inside getDataframe: df:', df.toString());
    return df;
  }

  public async init(): Promise<void> {
    const dataQueryList = (await grok.dapi.queries.list()).filter((dq) => dq.name.startsWith('Sequence'));
    const usersQuery = dataQueryList.find(
      (dq) => dq.name.endsWith('Users')
    )!;
    this.usersDf = await this.getDatafame(usersQuery);
    if (this.usersDf === undefined)
      throw new Error('users undefined!');
  };

  private usersDf?: DG.DataFrame;

  public async getUsers(): Promise<DG.DataFrame> {
    if (this.usersDf === undefined)
      await this.init();

    return this.usersDf!;
  }
}
