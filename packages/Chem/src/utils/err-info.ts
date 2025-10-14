
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export function errMsg(err: any): string {
  if (typeof err === 'string' || err instanceof String)
    return err as string;
  else if (err.constructor.name === 'StateError')
    return err['message'];
  else if (err.constructor.name === 'StateError' && '$thrownJsError' in err)
    return errMsg(err['$thrownJsError']);
  else if (err instanceof Error)
    return (err as Error).message;
  else
    return err.toString();
}

export function errStack(err: any): string | undefined {
  if (err instanceof Error)
    return err.stack;
  else if (err.constructor.name === 'StateError' && '$thrownJsError' in err)
    return errStack(err['$thrownJsError']);
  return undefined;
}

export function errInfo(err: any): [string, string | undefined] {
  return [errMsg(err), errStack(err)];
}
