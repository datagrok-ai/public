import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';


export function errorToConsole(err: any): string {
  if (typeof err === 'string' || err instanceof String) {
    return err as string;
  } else if ('$thrownJsError' in err) {
    return errorToConsole(err['$thrownJsError']);
  } else if (err instanceof Error) {
    return (err as Error).stack ?? (err as Error).message;
  } else {
    return err.toString();
  }
}

export function rectToConsole(rect: DG.Rect): string {
  return `(x=${rect.x}, y=${rect.y}, w=${rect.width}, h=${rect.height})`;
}
