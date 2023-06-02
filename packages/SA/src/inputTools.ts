// inputTools.ts

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export type VariedInputInfo = {
  name: string,
  type: DG.TYPE,
  min: any,
  max: any,
};