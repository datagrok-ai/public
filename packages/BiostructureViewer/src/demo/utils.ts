import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package';

export function handleError(err: any): void {
  const errMsg: string = err instanceof Error ? err.message : err.toString();
  const stack: string | undefined = err instanceof Error ? err.stack : undefined;
  grok.shell.error(errMsg);
  _package.logger.error(err.message, undefined, stack);
}
