/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';

import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {IMacromoleculeAssembler} from '@datagrok-libraries/chem-meta/src/macromolecule-assembler';
import {IMonomerLibBase} from '@datagrok-libraries/bio/src/types/monomer-library';
import {MolfileWithMap} from '@datagrok-libraries/bio/src/monomer-works/types';
import {getMolColName, getUnusedColName} from '@datagrok-libraries/bio/src/monomer-works/utils';
import {IHelmToMolfileConverter, ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {assembleRowWithMapSync, assembleSeqColumnSync} from '@datagrok-libraries/bio/src/monomer-works/assembler';
import {getMonomerLibKey} from '@datagrok-libraries/bio/src/monomer-works/library-json';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';

import {_package} from '../../../package';

/** A Macromolecule HELM column of the given HELM strings, for the assembler's input layer. */
export function helmListColumn(helmList: (string | null | undefined)[], name: string = 'helm'): DG.Column<string> {
  const col = DG.Column.fromStrings(name, Array.from(helmList, (helm) => helm ?? ''));
  col.semType = DG.SEMTYPE.MACROMOLECULE;
  col.meta.units = NOTATION.HELM;
  return col;
}

/** The column itself when it is a Macromolecule HELM column, else a HELM column of its values. */
export function asHelmColumn(helmCol: DG.Column<string>): DG.Column<string> {
  if (helmCol.semType === DG.SEMTYPE.MACROMOLECULE && helmCol.meta.units === NOTATION.HELM)
    return helmCol;
  return helmListColumn(Array.from(helmCol.toList()), helmCol.name);
}

/** HELM to molecule with the non-linear layout, through Chem's macromolecule assembler (Crux). The converter keeps
 * its synchronous interface (D14): it is made by the async SeqHelper.getHelmToMolfileConverter, which readies the
 * assembler's main-thread module and the library ({@link libKey}), and its methods run on that module with the
 * assembler's synchronous twins. Its own HELM parser, monomer wrappers and the RDKit and OpenChemLib steps are gone.
 * Stereo is always absolute: the assembler puts every specified centre in one ABS group. */
export class HelmToMolfileConverter implements IHelmToMolfileConverter {
  constructor(
    private readonly monomerLib: IMonomerLibBase,
    private readonly seqHelper: ISeqHelper,
    private readonly assembler: IMacromoleculeAssembler,
    /** The key {@link monomerLib} was readied with: the library as it was when the converter was made */
    private readonly libKey: string,
  ) { }

  /** The non-linear assembly of the column (the assembler lays the molecules out and keeps stereo absolute).
   * @param {DG.Column<string>} helmCol HELM
   * @param {boolean} _chiralityEngine no effect: stereo is always absolute
   * @param {RDModule} _rdKitModule deprecated and unused
   * @param {IMonomerLibBase} monomerLib the library to resolve with (the converter's when not given; another
   *   library must have been readied on the assembler too)
   * @return {DG.Column<string>} `molfile(<helmCol>)`, not added to the table
   * @deprecated Use SeqHelper.helmToAtomicLevel */
  convertToRdKitBeautifiedMolfileColumn(
    helmCol: DG.Column<string>, _chiralityEngine: boolean, _rdKitModule: RDModule | null, monomerLib: IMonomerLibBase,
  ): DG.Column<string> {
    const lib = monomerLib ?? this.monomerLib;
    const key = lib === this.monomerLib ? this.libKey : getMonomerLibKey(lib);
    const molCol = this.assembleColumn(helmCol, lib, key);
    molCol.name = getMolColName(helmCol.dataFrame, helmCol.name);
    return molCol;
  }

  public convertToMolfileV3KColumn(helmCol: DG.Column<string>): DG.Column<string> {
    const molCol = this.assembleColumn(helmCol, this.monomerLib, this.libKey);
    molCol.name = getUnusedColName(helmCol.dataFrame, `molfileV3K(${helmCol.name})`);
    return molCol;
  }

  /** Molfiles with their monomer maps (keyed by the HELM position, gaps counted), one row at a time on the
   * main-thread module; an empty one for an empty or failed HELM. */
  public convertToMolfileV3K(helmList: string[]): MolfileWithMap[] {
    const helmCol = helmListColumn(helmList);
    const resList: MolfileWithMap[] = new Array<MolfileWithMap>(helmList.length);
    for (let i = 0; i < helmList.length; ++i) {
      const helm = helmList[i];
      if (!helm) {
        resList[i] = MolfileWithMap.createEmpty();
        continue;
      }

      let resMolfile: MolfileWithMap;
      try {
        const res = assembleRowWithMapSync(helmCol, i, this.monomerLib, this.seqHelper, 'nonlinear',
          this.assembler, this.libKey);
        if (res.error)
          throw new Error(`Cannot get molfile of HELM '${helm}': ${res.error}`);
        resMolfile = new MolfileWithMap(res.molfile, res.monomers);
      } catch (err: any) {
        const [errMsg, errStack] = errInfo(err);
        _package.logger.error(errMsg, undefined, errStack);
        resMolfile = MolfileWithMap.createEmpty();
      }
      resList[i] = resMolfile;
    }
    return resList;
  }

  /** The non-linear molecule column of a HELM column, on the main-thread module (the legacy synchronous path). */
  private assembleColumn(helmCol: DG.Column<string>, lib: IMonomerLibBase, key: string): DG.Column<string> {
    const res = assembleSeqColumnSync(asHelmColumn(helmCol), lib, this.seqHelper, 'nonlinear', this.assembler, key);
    return res.molCol!;
  }
}
