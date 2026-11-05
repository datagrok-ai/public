import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {IMonomerLibBase} from '../types/monomer-library';

import {ISeqHandler} from './macromolecule/seq-handler';
import {MolfileWithMap} from '../monomer-works/types';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';

export type ToAtomicLevelRes = {
  molCol: DG.Column<string> | null,
  warnings: string[],
}

/** HELM to molecule (non-linear layout), served by Chem's macromolecule assembler. The methods are synchronous: the
 * converter comes from the async {@link ISeqHelper.getHelmToMolfileConverter}, which readies the assembler's
 * main-thread module and the library, and they run there. Stereo is always absolute (every specified centre in one
 * ABS group). */
export interface IHelmToMolfileConverter {
  /** @param rdKitModule deprecated and unused: the assembler lays the molecules out */
  convertToRdKitBeautifiedMolfileColumn(
    helmCol: DG.Column<string>, chiralityEngine: boolean, rdKitModule: RDModule | null, monomerLib: IMonomerLibBase
  ): DG.Column<string>;
  convertToMolfileV3KColumn(helmCol: DG.Column<string>): DG.Column<string>;
  /** Molfiles with their monomer maps, one per HELM (an empty one for an empty or failed HELM) */
  convertToMolfileV3K(helmList: string[]): MolfileWithMap[];
}

export interface ISeqHelper {
  /** HELM column to molecules with the non-linear layout (Chem's macromolecule assembler).
   * @param helmCol {DG.Column<string>} Macromolecules in Helm format
   * @param chiralityEngine {boolean} [true] No effect: stereo is always absolute
   * @param highlight {boolean} [false] Sets the result's `.sequence-src-highlight-monomers` tag
   * @param overrideMonomerLib {IMonomerLibBase} [null] Override monomer library for monomers from reactions
   **/
  helmToAtomicLevel(
    helmCol: DG.Column<string>, chiralityEngine?: boolean, highlight?: boolean, overrideMonomerLib?: IMonomerLibBase
  ): Promise<ToAtomicLevelRes>;

  /** One HELM to a molecule with its monomer map (non-linear layout), synchronously through a converter of
   * {@link getHelmToMolfileConverter}. chiralityEngine and beautifyMol have no effect: stereo is always absolute, the
   * layout is the depiction. */
  //sync method
  helmToAtomicLevelSingle(
    helm: string, converter: IHelmToMolfileConverter, chiralityEngine?: boolean, beautifyMol?: boolean
  ): MolfileWithMap;

  getSeqHandler(seqCol: DG.Column<string>): ISeqHandler;

  getSeqMonomers(seqCol: DG.Column<string>): string[];

  setUnitsToFastaColumn(sh: ISeqHandler): void;
  setUnitsToSeparatorColumn(sh: ISeqHandler): void;
  setUnitsToHelmColumn(sh: ISeqHandler): void;
  /** A converter whose synchronous methods use {@link monomerLib} (as it is now): readies the assembler first */
  getHelmToMolfileConverter(monomerLib: IMonomerLibBase): Promise<IHelmToMolfileConverter>
}

export async function getSeqHelper(): Promise<ISeqHelper> {
  const packageName = 'Bio';
  const funcList = DG.Func.find({package: packageName, name: `getSeqHelper`});
  if (funcList.length === 0)
    throw new Error(`Package '${packageName}' must be installed for SeqHelper.`);
  const res = (await funcList[0].prepare().call()).getOutputParamValue() as ISeqHelper;
  return res;
}
