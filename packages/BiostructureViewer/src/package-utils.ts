import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {IPdbHelper} from '@datagrok-libraries/bio/src/pdb/pdb-helper';
import {NglGlServiceBase} from '@datagrok-libraries/bio/src/viewers/ngl-gl-service';
import {NglGlDocService} from './utils/ngl-gl-doc-service';
import {PdbHelper} from './utils/pdb-helper';

// -- _package --

export class Package extends DG.Package {
  private _pLogger: DG.PackageLogger;

  get logger(): DG.PackageLogger {
    if (!this._pLogger) {
      this._pLogger = new class extends DG.PackageLogger {
        private logPrefix: string = 'BsV: ';

        constructor(_package: DG.Package) { super(_package); }

        debug(message: string, params?: object): void { super.debug(this.logPrefix + message, params); }
      }(this);
    }

    return this._pLogger;
  }
}

// -- _getNglGlService, _getPdbHelper--

type BsvWindowType = Window & {
  $pdbHelper?: IPdbHelper,
  $nglGlService?: NglGlServiceBase,
};
declare const window: BsvWindowType;

export function _getNglGlService() {
  if (!(window.$nglGlService)) {
    const svc: NglGlServiceBase = new NglGlDocService();
    window.$nglGlService = svc;
  }

  return window.$nglGlService;
}

export async function _getPdbHelper(): Promise<IPdbHelper> {
  if (!(window.$pdbHelper)) {
    const ph = await PdbHelper.getInstance(); // getPdbHelper
    window.$pdbHelper = ph;
  }

  return window.$pdbHelper;
}
