import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {ISeqHelper, ToAtomicLevelRes} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {getHelmHelper, IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {MolfileWithMap} from '@datagrok-libraries/bio/src/monomer-works/types';
import {getMolColName} from '@datagrok-libraries/bio/src/monomer-works/utils';
import {ChemTags} from '@datagrok-libraries/chem-meta/src/consts';
import {IMonomerLibBase} from '@datagrok-libraries/bio/src/types/index';

import {HelmToMolfileConverter} from '../helm-to-molfile/converter';
import {ISeqHandler} from '@datagrok-libraries/bio/src/utils/macromolecule/seq-handler';
import {SeqHandler} from './seq-handler';
import {Column} from 'datagrok-api/dg';
import {NOTATION, TAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {_package} from '../../package';

type SeqHelperWindowType = Window & { $seqHelperPromise?: Promise<SeqHelper> };
declare const window: SeqHelperWindowType;

export class SeqHelper implements ISeqHelper {
  constructor(
    private readonly libHelper: IMonomerLibHelper,
    private readonly rdKitModule: RDModule
  ) {}

  getSeqHandler(seqCol: DG.Column<string>): ISeqHandler {
    return SeqHandler.forColumn(seqCol, this);
  }

  getSeqMonomers(seqCol: Column<string>): string[] {
    const sh = this.getSeqHandler(seqCol);
    return Object.keys(sh.stats.freq);
  }

  // TODO: Move to the Helm package
  async getHelmToMolfileConverter(monomerLib: IMonomerLibBase): Promise<HelmToMolfileConverter> {
    const helmHelper: IHelmHelper = await getHelmHelper();
    return new HelmToMolfileConverter(helmHelper, this.rdKitModule, monomerLib);
  }

  helmToAtomicLevelSingle(
    helm: string, converter: HelmToMolfileConverter, chiralityEngine?: boolean, beautifyMol: boolean = true) {
    if (!helm)
      return MolfileWithMap.createEmpty();
    const molfileV3k = converter.convertToMolfileV3K([helm])[0];
    if (!molfileV3k || !molfileV3k.molfile)
      return MolfileWithMap.createEmpty();
    let mol: RDMol | null = null;
    try {
      let v3k = molfileV3k.molfile;
      if (beautifyMol) {
        mol = this.rdKitModule.get_mol(v3k);
        if (!mol)
          return MolfileWithMap.createEmpty();
        mol.set_new_coords();
        mol.normalize_depiction(1);
        mol.straighten_depiction(true);
        v3k = mol.get_v3Kmolblock();
      }
      if (chiralityEngine)
        v3k = converter.molV3KtoMolV3KOCL(v3k);
      return new MolfileWithMap(v3k, molfileV3k.monomers);
    } catch (err) {
      const [errMsg, errStack] = errInfo(err);
      _package.logger.error(errMsg, undefined, errStack);
      return MolfileWithMap.createEmpty();
    } finally {
      mol?.delete();
    }
  }

  async helmToAtomicLevel(
    helmCol: DG.Column<string>, chiralityEngine?: boolean, highlight?: boolean, overrideMonomerLib?: IMonomerLibBase
  ): Promise<ToAtomicLevelRes> {
    const monomerLib: IMonomerLibBase = overrideMonomerLib ?? this.libHelper.getMonomerLib();

    const df: DG.DataFrame = helmCol.dataFrame;
    const molColName: string = getMolColName(df, helmCol.name);

    const converter = await this.getHelmToMolfileConverter(monomerLib);

    //#endregion From HelmToMolfileConverter
    const helmList = helmCol.toList();
    const molList = new Array<string>(helmCol.length);
    const beautifyMolsChemFunc = DG.Func.find({package: 'Chem', name: 'beautifyMols'})[0];
    // both options set to false, ugly molfiles
    for (let i = 0; i < helmCol.length; i++) {
      molList[i] = (this.helmToAtomicLevelSingle(helmList[i], converter,
        chiralityEngine, !beautifyMolsChemFunc)).molfile;
    }
    // need to beautify the molfiles
    if (beautifyMolsChemFunc) {
      const beautifiedMols = await beautifyMolsChemFunc.apply({mols: molList});
      if (beautifiedMols && Array.isArray(beautifiedMols) && beautifiedMols.length === helmCol.length) {
        for (let i = 0; i < helmCol.length; i++)
          beautifiedMols[i] && (molList[i] = beautifiedMols[i]);
      }
    }
    //const molHlList = molfilesV3K.map((item: MolfileWithMap) => getMolHighlight(item.monomers.values(), monomerLib));

    const molCol = DG.Column.fromStrings(molColName, molList);
    molCol.semType = DG.SEMTYPE.MOLECULE;
    molCol.meta.units = DG.UNITS.Molecule.MOLBLOCK;
    molCol.setTag(ChemTags.SEQUENCE_SRC_COL, helmCol.name);

    return {molCol: molCol, warnings: []};
  }

  public setUnitsToFastaColumn(uh: SeqHandler) {
    if (uh.column.semType !== DG.SEMTYPE.MACROMOLECULE || uh.column.meta.units !== NOTATION.FASTA)
      throw new Error(`The column of notation '${NOTATION.FASTA}' must be '${DG.SEMTYPE.MACROMOLECULE}'.`);

    uh.column.meta.units = NOTATION.FASTA;
    SeqHandler.setTags(uh);
  }

  public setUnitsToSeparatorColumn(uh: SeqHandler, separator?: string) {
    if (uh.column.semType !== DG.SEMTYPE.MACROMOLECULE || uh.column.meta.units !== NOTATION.SEPARATOR)
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
