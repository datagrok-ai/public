import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {v4 as uuidv4} from 'uuid';

import {errInfo} from './err-info';

export enum LogLevel { error = 0, warning = 1, info = 2, debug = 3}

export interface ILogger {
  error(message: any, params?: object | undefined, stackTrace?: string | undefined): void;
  warning(message: string, params?: object | undefined): void;
  info(message: string, params?: object | undefined): void;
  debug(message: string, params?: object | undefined): void;
}

export class PromiseSyncer {
  private promise: Promise<void> = Promise.resolve();

  constructor(
    private readonly logger: ILogger,
  ) {}

  private syncCounter: number = 0;

  public sync(logPrefix: string, func: () => Promise<void>): void {
    const syncId: number = ++this.syncCounter;
    this.logger.debug(`${logPrefix}, SYNC syncId = ${syncId}, IN `);
    this.promise = this.promise.then(async () => {
      this.logger.debug(`${logPrefix}, SYNC syncId = ${syncId}, START `);
      await func();
      this.logger.debug(`${logPrefix}, SYNC syncId = ${syncId}, END `);
    }).catch((err: any) => {
      const [errMsg, errStack] = errInfo(err);
      this.logger.error(`${logPrefix}, SYNC syncId = ${syncId}, ERROR:\n${errMsg}`, undefined, errStack);
    });
    // this.logger.debug(`${logPrefix}, SYNC syncId = ${syncId}, OUT`);
  }
}
