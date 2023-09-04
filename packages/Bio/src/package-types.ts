import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {Observable, Subject} from 'rxjs';
import {ObjectPropertyBag} from 'datagrok-api/dg';

/** Names of package properties/settings declared in properties section of {@link './package.json'} */
export const enum BioPackagePropertiesNames {
  MaxMonomerLength = 'MaxMonomerLength',
  TooltipWebLogo = 'TooltipWebLogo',
  DefaultSeparator = 'DefaultSeparator',
}


export class BioPackageProperties extends Map<string, any> {
  private _onPropertyChanged: Subject<string> = new Subject<string>();
  public get onPropertyChanged(): Observable<string> { return this._onPropertyChanged; }

  /** Monomer name maximum length displayed in short mode. */
  public get MaxMonomerLength(): number {
    return super.get(BioPackagePropertiesNames.MaxMonomerLength) as number;
  }

  public set MaxMonomerLength(value: number) {
    super.set(BioPackagePropertiesNames.MaxMonomerLength, value);
    this._onPropertyChanged.next(BioPackagePropertiesNames.MaxMonomerLength);
  }

  public get TooltipWebLogo(): boolean {
    return super.get(BioPackagePropertiesNames.TooltipWebLogo) as boolean;
  }

  public set TooltipWebLogo(value: boolean) {
    super.set(BioPackagePropertiesNames.TooltipWebLogo, value);
    this._onPropertyChanged.next(BioPackagePropertiesNames.TooltipWebLogo);
  }

  public get DefaultSeparator(): string {
    return super.get(BioPackagePropertiesNames.DefaultSeparator) as string;
  }

  public set DefaultSeparator(value: string) {
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
  /** Package properties/settings declared in properties section of {@link './package.json'} */
  public get properties(): BioPackageProperties { return this._properties; };

  public set properties(value: BioPackageProperties) { this._properties = value; }
}
