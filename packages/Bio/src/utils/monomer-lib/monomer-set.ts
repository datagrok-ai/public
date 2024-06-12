import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  IMonomerLib, IMonomerLink, IMonomerLinkData, IMonomerSet, IMonomerSetPlaceholder, Monomer, MonomerType, PolymerType
} from '@datagrok-libraries/bio/src/types';

export class MonomerSetPlaceholder implements IMonomerSetPlaceholder {
  public readonly monomers: Monomer[];

  constructor(
    private readonly monomerLib: IMonomerLib,
    public symbol: string,
    public readonly polymerType: PolymerType,
    public readonly monomerType: MonomerType,
    public readonly monomerLinks: IMonomerLinkData[],
  ) {
    this.monomers = this.monomerLinks.map((mLink) => {
      const resM = this.monomerLib.getMonomer(this.polymerType, mLink.symbol);
      if (!resM)
        throw new Error('Monomer not found: ');
      if (resM.lib?.source != mLink.libName)
        throw new Error('Monomer found in different library.');
      return resM;
    });
  }
}

export class MonomerSet implements IMonomerSet {
  public constructor(
    public readonly description: string,
    public readonly placeholders: IMonomerSetPlaceholder[],
    public readonly source: string | undefined = undefined,
    public readonly error: string | undefined = undefined,
  ) {}
}
