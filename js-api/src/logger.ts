import {LOG_LEVEL} from './const';
import {toDart} from './wrappers';
import {Package} from './entities';
import {IDartApi} from "./api/grok_api.g";
import dayjs from "dayjs";
import {Accordion} from "./widgets";
import {toJs} from "./wrappers";
import {DataFrame} from "./dataframe";

const api: IDartApi = <any>window;

export type LogMessage = {level: LOG_LEVEL, message: string, params?: object, type?: string, stackTrace?: string};

type LoggerPutCallback = (logRecord: LogMessage) => void;


export class Logger {
  putCallback?: LoggerPutCallback;
  private readonly dart: any;
  private static consoleLogs: object[] = [];

  constructor(putCallback?: LoggerPutCallback, options?: {staticLogger?: boolean, params?: object, dartLogger?: any}) {
    this.putCallback = putCallback;
    if (options?.staticLogger != true)
      this.dart = options?.dartLogger ?? api.grok_GetLogger(toDart(options?.params));
  }

  static create(options?: {params: object}) {
    return new Logger(undefined, {params: options?.params});
  }

  static getStatic() {
    return new Logger(undefined, {staticLogger: true});
  }

  static translateStackTrace(stackTrace: string): Promise<string> {
     return api.grok_Log_TranslateStackTrace(stackTrace);
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
  error(message: any, params?: object, stackTrace?: string): void {
    if (message instanceof Error) {
      stackTrace ??= message.stack;
      message = message.message;
    }
    this._log({level: LOG_LEVEL.ERROR, message, params, stackTrace});
  }

  _log(msg: LogMessage) {
    if (this.putCallback != null)
      this.putCallback(msg);
    msg.stackTrace ??= new Error().stack;
    api.grok_Log(this.dart, msg.level, msg.message, toDart(msg.params), msg.type, msg.stackTrace);
  }

  static getConsoleOutput(): string {
    return JSON.stringify(this.consoleLogs, null, 2);
  }

  static interceptConsoleOutput(): void {
    const regex = new RegExp('(?:\\d{4}-\\d{2}-\\d{2})?(?:[ T]\\d{2}:\\d{2}:\\d{2})?[.,]\\d{3}');
    const intercept = (f: (message?: any, ...optionalParams: any[]) => void) => {
      const std = f.bind(console);
      return (...args: any[]) => {
        try {
          let message = args.map((x) => `${x}`).join(' ');
          let time = dayjs().utc().valueOf();
          if (regex.test(message)) {
            const result = regex.exec(message);
            if (result && result.length > 0) {
              time = dayjs(result[0]).utc().valueOf();
              message = message.replace(regex, '');
            }
          }
          this.consoleLogs.push({'time': time, 'message': message});
        } catch (_) {}
        std(...args);
      }
    };

    console.log = intercept(console.log);
    console.info = intercept(console.info);
    console.warn = intercept(console.warn);
    console.error = intercept(console.error);
  }

  static set autoReportOptions(options: {[key: string]: any}) {
    api.grok_Set_AutoReport_Options(toDart(options));
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
    //@ts-ignore
    msg.params['packageName'] = this.package.name;
    super._log(msg);
  }
}

