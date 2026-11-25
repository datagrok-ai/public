import * as grok from 'datagrok-api/grok';
import {errorMessage, ForgeError} from '../forge-error';
import {_package} from '../package';

export function reportError(e: unknown, isRepeating: boolean = false): void {
  const message = errorMessage(e);
  const options = isRepeating ? {oneTimeKey: message} : undefined;
  if (e instanceof ForgeError)
    grok.shell.warning(message, options);
  else {
    grok.shell.error(message, options);
    _package.logger.error(e);
  }
}
