import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import type {MonomerType, PolymerType} from '@datagrok-libraries/bio/src/helm/types';
import {
  IMonomerLib, IMonomerLinkData, IMonomerSet, IMonomerSetPlaceholder, Monomer
} from '@datagrok-libraries/bio/src/types';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';

export class MonomerSetPlaceholder implements IMonomerSetPlaceholder {
  public readonly monomers: Monomer[];

  public readonly error: string | null = null;

  constructor(
    private readonly monomerLib: IMonomerLib,
    public symbol: string,
    public readonly polymerType: PolymerType,
    public readonly monomerType: MonomerType,
    public readonly monomerLinks: IMonomerLinkData[],
  ) {
    try {
      this.monomers = this.monomerLinks.map((mLink) => {
        const resM = this.monomerLib.getMonomer(this.polymerType, mLink.symbol);
        if (!resM)
          throw new Error('Monomer not found: ');
        if (resM.lib?.source != mLink.source)
          throw new Error(`Monomer '${symbol}' found in different library.`);
        return resM;
      });
    } catch (err: any) {
      const [errMsg, errStack] = errInfo(err);
      this.error = errMsg;
      this.monomers = [];
    }
  }
}

export class MonomerSet implements IMonomerSet {
  public constructor(
    public readonly description: string,
    public placeholders: IMonomerSetPlaceholder[],
    public readonly source: string | undefined = undefined,
    public readonly error: string | undefined = undefined,
  ) {}

  updateSets(setList: IMonomerSet[], reload: boolean = false): void {
    if (reload)
      this.placeholders = [];
    for (const _set of setList)
      if (!_set.error) this._updateSetInt(_set);

    // TODO: File onChanged
  }

  private _updateSetInt(set: IMonomerSet): void {
    for (const setPh of set.placeholders)
      this.placeholders.push(setPh);
  }
}
