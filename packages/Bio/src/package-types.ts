import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {Observable, Subject} from 'rxjs';

import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {LoggerWrapper} from '@datagrok-libraries/bio/src/utils/logger';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';

/** Names of package properties/settings declared in properties section of {@link './package.json'} */
export const enum BioPackagePropertiesNames {
  MonomerWidthMode = 'MonomerWidthMode',
  MaxMonomerLength = 'MaxMonomerLength',
  TooltipWebLogo = 'TooltipWebLogo',
  DefaultSeparator = 'DefaultSeparator',
}


export class BioPackageProperties extends Map<string, any> {
  private _onPropertyChanged: Subject<string> = new Subject<string>();
  public get onPropertyChanged(): Observable<string> { return this._onPropertyChanged; }

  /** Monomer symbol maximum length displayed, null for unlimited. */
  public get maxMonomerLength(): number | null {
    const vs = super.get(BioPackagePropertiesNames.MaxMonomerLength);
    return vs === 'long' ? null : parseInt(vs);
  }

  public set maxMonomerLength(value: number | null) {
    const vs = value === null ? 'long' : value.toString();
    super.set(BioPackagePropertiesNames.MaxMonomerLength, vs);
    this._onPropertyChanged.next(BioPackagePropertiesNames.MaxMonomerLength);
  }

  public get tooltipWebLogo(): boolean {
    return super.get(BioPackagePropertiesNames.TooltipWebLogo) as boolean;
  }

  public set tooltipWebLogo(value: boolean) {
    super.set(BioPackagePropertiesNames.TooltipWebLogo, value);
    this._onPropertyChanged.next(BioPackagePropertiesNames.TooltipWebLogo);
  }

  public get defaultSeparator(): string {
    return super.get(BioPackagePropertiesNames.DefaultSeparator) as string;
  }

  public set defaultSeparator(value: string) {
    if (value.length !== 1) throw new Error('The separator must be of length one.');
    super.set(BioPackagePropertiesNames.DefaultSeparator, value);
    this._onPropertyChanged.next(BioPackagePropertiesNames.DefaultSeparator);
  }

  constructor(source: any) {
    super(Object.entries(source));
  }
}

export class BioPackage extends DG.Package {
  private _properties: BioPackageProperties;

  private _rdKitModule: RDModule;
  public get rdKitModule(): RDModule { return this._rdKitModule;};

  /** Package properties/settings declared in properties section of {@link './package.json'} */
  public get properties(): BioPackageProperties { return this._properties; };

  public set properties(value: BioPackageProperties) { this._properties = value; }

  private _initialized: boolean = false;

  public get initialized(): boolean { return this._initialized; }

  constructor(opts: { debug: boolean } = {debug: false}) {
    super();
    // @ts-ignore
    super._logger = new LoggerWrapper(super.logger, opts.debug);
  }

  public completeInit(rdKitModule: RDModule): void {
    this._rdKitModule = rdKitModule;
    this._initialized = true;
  }

  handleErrorUI(err: any) {
    const [errMsg, errStack] = errInfo(err);
    grok.shell.error(errMsg);
    this.logger.error(errMsg, undefined, errStack);
  }
}
