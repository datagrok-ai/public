/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package';

export abstract class DataLoaderBase {
  abstract get users(): DG.DataFrame;
  // abstract get users(): DG.DataFrame;
  // abstract get users(): DG.DataFrame;
  // abstract get users(): DG.DataFrame;
  // abstract get users(): DG.DataFrame;
}

export class DataLoaderDB extends DataLoaderBase {
  public async init(): Promise<void> {
    const dataQueryList = (await grok.dapi.queries.list()).filter((dq) => dq.name.startsWith('Sequence'));
    const usersQuery = dataQueryList.find(
      (dq) => dq.name.endsWith('Users')
    )!;
    [this.usersDf] = await Promise.all([
      (async () => {
        const c = await usersQuery.prepare().call();
        const df: DG.DataFrame = c.getOutputParamValue();
        return df;
      })(), // guaranteed to be executed
      // todo: 4 more for othre queries
    ]);
  };

  private usersDf?: DG.DataFrame;

  get users(): DG.DataFrame { return this.usersDf!; }
}
