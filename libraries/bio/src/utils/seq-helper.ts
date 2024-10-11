import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {IMonomerLibBase} from '../types';

import {ISeqHandler} from './macromolecule/seq-handler';

export type ToAtomicLevelRes = {
  molCol: DG.Column<string> | null,
  warnings: string[],
}

export interface ISeqHelper {
  /**
   * @param helmCol {DG.Column<string>} Macromolecules in Helm format
   * @param chiralityEngine {boolean} [true] Use chirality engine for molecule visualization
   * @param highlight {boolean} [true] Generates molHighlightCol of result
   * @param overrideMonomerLib {IMonomerLibBase} [null] Override monomer library for monomers from reactions
   **/
  helmToAtomicLevel(
    helmCol: DG.Column<string>, chiralityEngine?: boolean, highlight?: boolean, overrideMonomerLib?: IMonomerLibBase
  ): Promise<ToAtomicLevelRes>;

  getSeqHandler(seqCol: DG.Column<string>): ISeqHandler;

  getSeqMonomers(seqCol: DG.Column<string>): string[];
}

export async function getSeqHelper(): Promise<ISeqHelper> {
  const packageName = 'Bio';
  const funcList = DG.Func.find({package: packageName, name: `getSeqHelper`});
  if (funcList.length === 0)
    throw new Error(`Package '${packageName}' must be installed for SeqHelper.`);
  const res = (await funcList[0].prepare().call()).getOutputParamValue() as ISeqHelper;
  return res;
}
