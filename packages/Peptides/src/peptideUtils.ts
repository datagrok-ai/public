import {getSeqHelper, ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';

export class PeptideUtils {
  private static _secHelper: ISeqHelper;
  public static getSeqHelper(): ISeqHelper {
    if (!this._secHelper)
      throw new Error('SeqHelper is not initialized');
    return this._secHelper;
  }

  public static async loadSeqHelper(): Promise<void> {
    this._secHelper ??= await getSeqHelper();
  }
}
