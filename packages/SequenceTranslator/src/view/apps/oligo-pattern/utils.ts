/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import { STRANDS, StrandType } from '../../../model/pattern-app/const';

// import {axolabsStyleMap} from '../../../model/data-loading-utils/json-loader';

export function applyToAllStrands(callback: (strand: StrandType) => any) {
  const result = Object.fromEntries(STRANDS.map((strand) => callback(strand)));
  return result;
}
