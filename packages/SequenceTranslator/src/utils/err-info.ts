import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';

import {_package} from '../package';

export function defaultErrorHandler(err: any, shell: boolean = true): [string, string | undefined] {
  const [errMsg, errStack] = errInfo(err);
  _package.logger.error(errMsg, undefined, errStack);
  if (shell) grok.shell.error(errMsg);
  return [errMsg, errStack];
}
