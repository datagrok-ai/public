import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export interface ILogger {
  error(message: any, params?: object | undefined, stackTrace?: string | undefined): void;
  warning(message: string, params?: object | undefined): void;
  info(message: string, params?: object | undefined): void;
  debug(message: string, params?: object | undefined): void;
}
