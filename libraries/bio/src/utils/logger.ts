import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export interface ILogger {
  error(message: any, params?: object | undefined, stackTrace?: string | undefined): void;
  warning(message: string, params?: object | undefined): void;
  info(message: string, params?: object | undefined): void;
  debug(message: string, params?: object | undefined): void;
}

export class LoggerWrapper implements ILogger {
  constructor(
    private readonly target: ILogger,
    private debugEnabled: boolean
  ) {}

  error(message: any, params?: object | undefined, stackTrace?: string | undefined): void {
    return this.target.error(message, params, stackTrace);
  }

  warning(message: string, params?: object | undefined): void {
    return this.target.warning(message, params);
  }

  info(message: string, params?: object | undefined): void {
    return this.target.info(message, params);
  }

  debug(message: string, params?: object | undefined): void {
    if (this.debugEnabled)
      return this.target.debug(message, params);
  }
}
