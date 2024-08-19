import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

export type DataProviderFunc<TId, TValue> = (id: TId) => Promise<TValue>;

export async function getDataProviderList(semType: string): Promise<DG.Func[]> {
  return grok.dapi.functions
    .include('params,entityTags')
    .filter(`options.dataProvider="${semType}"`)
    .list();
}
