import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';

// -- _package --

export class DockingPackage extends DG.Package {
  private _pLogger!: DG.PackageLogger;

  get logger(): DG.PackageLogger {
    if (!this._pLogger) {
      this._pLogger = new class extends DG.PackageLogger {
        private logPrefix: string = 'Docking: ';

        constructor(_package: DG.Package) { super(_package); }

        debug(message: string, params?: object): void { super.debug(this.logPrefix + message, params); }
      }(this);
    }

    return this._pLogger;
  }

  handleErrorUI(err: any) {
    const [errMsg, errStack] = errInfo(err);
    grok.shell.error(errMsg);
    this.logger.error(errMsg, undefined, errStack);
  }
}
