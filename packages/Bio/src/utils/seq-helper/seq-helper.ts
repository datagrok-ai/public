import * as DG from 'datagrok-api/dg';

import {ISeqHelper, ToAtomicLevelRes} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {MolfileWithMap} from '@datagrok-libraries/bio/src/monomer-works/types';
import {getMolColName} from '@datagrok-libraries/bio/src/monomer-works/utils';
import {ChemTags} from '@datagrok-libraries/chem-meta/src/consts';
import {IMonomerLibBase, IMonomerLibHelper} from '@datagrok-libraries/bio/src/types/monomer-library';

import {
  assembleSeqColumn, getMacromoleculeAssembler, readyAssemblerLibrary,
} from '@datagrok-libraries/bio/src/monomer-works/assembler';

import {asHelmColumn, HelmToMolfileConverter} from '../helm-to-molfile/converter/converter';
import {ISeqHandler} from '@datagrok-libraries/bio/src/utils/macromolecule/seq-handler';
import {SeqHandler} from './seq-handler';
import {Column} from 'datagrok-api/dg';
import {NOTATION, TAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';

export class SeqHelper implements ISeqHelper {
  constructor(private readonly libHelper: IMonomerLibHelper) {}

  getSeqHandler(seqCol: DG.Column<string>): ISeqHandler {
    return SeqHandler.forColumn(seqCol, this);
  }

  getSeqMonomers(seqCol: Column<string>): string[] {
    const sh = this.getSeqHandler(seqCol);
    return Object.keys(sh.stats.freq);
  }

  /** A converter for {@link monomerLib} as it is now: the assembler's main-thread module and the library are readied
   * here, so the converter's methods (and {@link helmToAtomicLevelSingle}) can be synchronous. */
  async getHelmToMolfileConverter(monomerLib: IMonomerLibBase): Promise<HelmToMolfileConverter> {
    const assembler = await getMacromoleculeAssembler();
    const libKey = await readyAssemblerLibrary(assembler, monomerLib);
    return new HelmToMolfileConverter(monomerLib, this, assembler, libKey);
  }

  /** One HELM to a molecule with its monomer map, non-linear layout, synchronously through a converter of
   * {@link getHelmToMolfileConverter}. chiralityEngine and beautifyMol have no effect: stereo is always absolute,
   * and the layout is the assembler's depiction. */
  helmToAtomicLevelSingle(
    helm: string, converter: HelmToMolfileConverter, _chiralityEngine?: boolean, _beautifyMol: boolean = true,
  ): MolfileWithMap {
    if (!helm)
      return MolfileWithMap.createEmpty();
    return converter.convertToMolfileV3K([helm])[0] ?? MolfileWithMap.createEmpty();
  }

  /** HELM column to molecules with the non-linear layout, in the assembler's workers. Rows that cannot be built
   * give '' and a warning each. chiralityEngine has no effect (stereo is always absolute); highlight sets the
   * `.sequence-src-highlight-monomers` tag of the result. */
  async helmToAtomicLevel(
    helmCol: DG.Column<string>, chiralityEngine?: boolean, highlight?: boolean, overrideMonomerLib?: IMonomerLibBase
  ): Promise<ToAtomicLevelRes> {
    const monomerLib: IMonomerLibBase = overrideMonomerLib ?? this.libHelper.getMonomerLib();
    const res = await assembleSeqColumn(asHelmColumn(helmCol), monomerLib, this, 'nonlinear');
    if (res.molCol) {
      res.molCol.name = getMolColName(helmCol.dataFrame, helmCol.name);
      res.molCol.setTag(ChemTags.SEQUENCE_SRC_COL, helmCol.name);
      res.molCol.setTag(ChemTags.SEQUENCE_SRC_HL_MONOMERS, String(!!highlight));
    }
    return res;
  }

  public setUnitsToFastaColumn(uh: SeqHandler) {
    if (uh.column.semType !== DG.SEMTYPE.MACROMOLECULE || uh.column.meta.units !== NOTATION.FASTA)
      throw new Error(`The column of notation '${NOTATION.FASTA}' must be '${DG.SEMTYPE.MACROMOLECULE}'.`);

    uh.column.meta.units = NOTATION.FASTA;
    SeqHandler.setTags(uh);
  }

  public setUnitsToSeparatorColumn(uh: SeqHandler, separator?: string) {
    if (uh.column.semType !== DG.SEMTYPE.MACROMOLECULE)
      throw new Error(`The column of notation '${NOTATION.SEPARATOR}' must be '${DG.SEMTYPE.MACROMOLECULE}'.`);
    if (!separator)
      throw new Error(`The column of notation '${NOTATION.SEPARATOR}' must have the separator tag.`);

    uh.column.meta.units = NOTATION.SEPARATOR;
    uh.column.setTag(TAGS.separator, separator);
    SeqHandler.setTags(uh);
  }

  public setUnitsToHelmColumn(uh: SeqHandler) {
    if (uh.column.semType !== DG.SEMTYPE.MACROMOLECULE)
      throw new Error(`The column of notation '${NOTATION.HELM}' must be '${DG.SEMTYPE.MACROMOLECULE}'`);

    uh.column.meta.units = NOTATION.HELM;
    SeqHandler.setTags(uh);
  }
}
