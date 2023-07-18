import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {Observable, Subject} from 'rxjs';
import {ObjectPropertyBag} from 'datagrok-api/dg';

/** Names of package properties/settings declared in properties section of {@link './package.json'} */
export const enum BioPackagePropertiesNames {
  MaxMonomerLength = 'MaxMonomerLength',
  TooltipWebLogo = 'TooltipWebLogo',
}


export class BioPackageProperties extends Map<string, any> {

  private _onPropertyChanged: Subject<string> = new Subject<string>();
  public get onPropertyChanged(): Observable<string> { return this._onPropertyChanged; }

  /** Monomer name maximum length displayed in short mode. */
  public get maxMonomerLength(): number {
    return super.get(BioPackagePropertiesNames.MaxMonomerLength) as unknown as number;
  }

  public set maxMonomerLength(value: number) {
    super.set(BioPackagePropertiesNames.MaxMonomerLength, value as unknown as object);
    this._onPropertyChanged.next(BioPackagePropertiesNames.MaxMonomerLength);
  }

  public get tooltipWebLogo(): boolean {
    return super.get(BioPackagePropertiesNames.TooltipWebLogo) as unknown as boolean;
  }

  public set tooltipWebLogo(value: boolean) {
    super.set(BioPackagePropertiesNames.TooltipWebLogo, value as unknown as boolean);
    this._onPropertyChanged.next(BioPackagePropertiesNames.TooltipWebLogo);
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
