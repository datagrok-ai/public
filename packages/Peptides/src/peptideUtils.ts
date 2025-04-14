import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types';
import {getSeqHelper, ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';

export class PeptideUtils {
  private static _secHelper: ISeqHelper;
  private static _monomerLibHelper: IMonomerLibHelper;
  public static getMonomerLib(): IMonomerLib {
    if (!this._monomerLibHelper)
      throw new Error('MonomerLib is not initialized');
    return this._monomerLibHelper.getMonomerLib();
  }
  public static getSeqHelper(): ISeqHelper {
    if (!this._secHelper)
      throw new Error('SeqHelper is not initialized');
    return this._secHelper;
  }

  public static async loadSeqHelper(): Promise<void> {
    this._secHelper ??= await getSeqHelper();
  }

  public static async loadComponents(): Promise<void> {
    await this.loadSeqHelper();
    await this.loadMonomerLib();
  }

  public static async loadMonomerLib(): Promise<void> {
    this._monomerLibHelper ??= await getMonomerLibHelper();
    await this._monomerLibHelper.awaitLoaded();
  }
}
