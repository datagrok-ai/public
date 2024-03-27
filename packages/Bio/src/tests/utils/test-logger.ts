import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';


export class TestLogger implements ILogger {
  errorList: { message: any, params?: object, stackTrace?: string }[] = [];
  warningList: { message: string, params?: object }[] = [];
  infoList: { message: string, params?: object }[] = [];
  debugList: { message: string, params?: object }[] = [];

  error(message: any, params?: object, stackTrace?: string): void {
    this.errorList.push({message, params, stackTrace});
  }

  warning(message: string, params?: object): void {
    this.warningList.push({message, params});
  }

  info(message: string, params?: object): void {
    this.infoList.push({message, params});
  }

  debug(message: string, params?: object): void {
    this.debugList.push({message, params});
  }
}
