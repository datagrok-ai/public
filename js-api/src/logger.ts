import {LOG_LEVEL} from './const';
import {toDart} from './wrappers';
import {Package} from './entities';

let api = <any>window;

export type LogMessage = {level: LOG_LEVEL, message: string, params?: object, type?: string, stackTrace?: string};

type LoggerPutCallback = (logRecord: LogMessage) => void;


export class Logger {
  putCallback?: LoggerPutCallback;
  private readonly dart: any;

  constructor(putCallback?: LoggerPutCallback, options?: {staticLogger?: boolean, params?: object}) {
    this.putCallback = putCallback;
    if (options?.staticLogger != true)
      this.dart = api.grok_GetLogger(toDart(options?.params));
  }

  static create(options?: {params: object}) {
    return new Logger(undefined, {params: options?.params});
  }

  static getStatic() {
    return new Logger(undefined, {staticLogger: true});
  }

  /** @Obsolete, for backward compatibility, use {audit} instead **/
  log(message: string, params: object, type: string = 'log'): void {
    this._log({level: LOG_LEVEL.AUDIT, message, params, type});
  }

  /** Reports audit record to Datagrok **/
  audit(message: string, params?: object, type: string = 'log'): void {
    this._log({level: LOG_LEVEL.AUDIT, message, params, type});
  }

  /** Reports usage record to Datagrok **/
  usage(message: string, params?: object, type: string = 'usage'): void {
    this._log({level: LOG_LEVEL.USAGE, message, params, type});
  }

  /** Reports info record to Datagrok **/
  info(message: string, params?: object): void {
    this._log({level: LOG_LEVEL.INFO, message, params});
  }

  /** Reports debug record to Datagrok **/
  debug(message: string, params?: object): void {
    this._log({level: LOG_LEVEL.DEBUG, message, params});
  }

  /** Reports warning record to Datagrok **/
  warning(message: string, params?: object): void {
    this._log({level: LOG_LEVEL.WARNING, message, params});
  }

  /** Reports error record to Datagrok **/
  error(message: string, params?: object, stackTrace?: string): void {
    this._log({level: LOG_LEVEL.ERROR, message, params, stackTrace});
  }

  _log(msg: LogMessage) {
    if (this.putCallback != null)
      this.putCallback(msg);
    api.grok_Log(this.dart, msg.level, msg.message, toDart(msg.params), msg.type, msg.stackTrace);
  }
}


export class PackageLogger extends Logger {
  private package: Package;

  constructor(_package: Package) {
    super();
    this.package = _package;
  }

  _log(msg: LogMessage) {
    msg.params ??= {};
    //@ts-ignore
    msg.params['package'] = this.package.dart;
    super._log(msg);
  }
}

