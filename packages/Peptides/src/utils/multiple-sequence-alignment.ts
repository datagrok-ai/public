/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
// eslint-disable-next-line no-unused-vars
import * as grok from 'datagrok-api/grok';
// eslint-disable-next-line no-unused-vars
import * as ui from 'datagrok-api/ui';
// eslint-disable-next-line no-unused-vars
import * as DG from 'datagrok-api/dg';

import biomsa from 'biomsa';

import {AlignedSequenceEncoder} from '@datagrok-libraries/utils/src/sequence-encoder';

export async function doMSA(col: DG.Column): Promise<DG.Column> {
  const sequences = col.toList().map((v: string, _) => AlignedSequenceEncoder.clean(v).replace('-', ''));
  const aligned = await biomsa.align(sequences);
  return DG.Column.fromStrings('aligned', aligned);
}
