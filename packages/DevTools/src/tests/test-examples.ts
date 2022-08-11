import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {list} from 'datagrok-api/ui';

export function throwsException(): void {
  throw 'An exception';
}

export function returnsFine(): void {
}
