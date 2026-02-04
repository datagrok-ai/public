import * as DG from 'datagrok-api/dg';

import {delay} from '@datagrok-libraries/test/src/test';


export async function getOptions(viewer: DG.Viewer): Promise<Indexable> {
  await delay(100);
  const options = (viewer.getOptions(true) as {id: string, type: string, look: Indexable}).look;
  delete options['#type'];
  delete options.table;
  return options;
}

export interface Indexable { [key: string]: any }
