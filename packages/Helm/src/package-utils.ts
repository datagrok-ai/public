/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';

import {Unsubscribable} from 'rxjs';

import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types/monomer-library';
import {LoggerWrapper} from '@datagrok-libraries/bio/src/utils/logger';
import {IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {IMonomerLibHelper} from '@datagrok-libraries/bio/src/types/monomer-library';
import {ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';

import * as DG from 'datagrok-api/dg';

import {_package} from './package';

// hwe migration (Phase 7): the legacy Dojo / JSDraw2 / HELMWebEditor loader and
// the `org.helm.webeditor.Monomers` patching that lived here are gone. The
// standalone `@datagrok-libraries/hwe` editor is self-contained (no global
// scripts to inject) and reads Bio's monomer library directly via
// `bridgeMonomerLib`, so package init just wires the helpers together.

export class HelmPackage extends DG.Package {
  private _helmHelper?: IHelmHelper;
  public get helmHelper(): IHelmHelper {
    if (!this._helmHelper)
      throw new Error('Package Helm .helmHelper is not initialized');
    return this._helmHelper;
  };

  public get seqHelper(): ISeqHelper {
    return this.helmHelper!.seqHelper;
  }

  public _libHelper!: IMonomerLibHelper;

  constructor(opts: { debug: boolean } = {debug: false}) {
    super();
    // @ts-ignore
    super._logger = new LoggerWrapper(super.logger, opts.debug);
  }

  // -- MonomerLib --

  public get monomerLib(): IMonomerLib {
    if (!this._libHelper)
      throw new Error(`Helm: _package.libHelper is not initialized yet`);
    return this._libHelper.getMonomerLib();
  }

  private _monomerLibSub?: Unsubscribable;

  private _initialized: boolean = false;

  /** Requires Bio initialized (monomer library). */
  completeInit(helmHelper: IHelmHelper, libHelper: IMonomerLibHelper): void {
    this._helmHelper = helmHelper;
    this._libHelper = libHelper;

    const lib = this.monomerLib;
    this._monomerLibSub = lib.onChanged
      .subscribe(this.monomerLibOnChangedHandler.bind(this));

    this._initialized = true;
  }

  monomerLibOnChangedHandler(): void {
    const logPrefix = `Helm: _package.monomerLibOnChangedHandler()`;
    try {
      const libSummary = this.monomerLib!.getSummaryObj();
      const isLibEmpty = Object.keys(libSummary).length == 0;
      const libSummaryLog = isLibEmpty ? 'empty' : Object.entries(libSummary)
        .map(([pt, count]) => `${pt}: ${count}`)
        .join(', ');
      _package.logger.debug(`${logPrefix}, start, lib: { ${libSummaryLog} }`);

      const libSummaryHtml = isLibEmpty ? 'empty' : Object.entries(libSummary)
        .map(([pt, count]) => `${pt} ${count}`)
        .join('<br />');
      const libMsg: string = `Monomer lib updated:<br /> ${libSummaryHtml}`;
      grok.shell.info(libMsg);
    } catch (err: any) {
      const errMsg = errorToConsole(err);
      console.error(`${logPrefix} error:\n` + errMsg);
      // throw err; // Prevent disabling event handler
    }
  }
}
