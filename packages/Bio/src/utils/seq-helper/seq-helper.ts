import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {ISeqHelper, ToAtomicLevelRes} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {getHelmHelper, IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {MolfileWithMap} from '@datagrok-libraries/bio/src/monomer-works/types';
import {getMolColName, hexToPercentRgb} from '@datagrok-libraries/bio/src/monomer-works/utils';
import {ChemTags} from '@datagrok-libraries/chem-meta/src/consts';
import {getMolHighlight} from '@datagrok-libraries/bio/src/monomer-works/seq-to-molfile';
import {IMonomerLibBase} from '@datagrok-libraries/bio/src/types/index';

import {HelmToMolfileConverter} from '../helm-to-molfile/converter';
import {MonomerLibManager} from '../monomer-lib/lib-manager';

import {_package, getMonomerLibHelper} from '../../package';

type SeqHelperWindowType = Window & { $seqHelperPromise?: Promise<SeqHelper> };
declare const window: SeqHelperWindowType;

export class SeqHelper implements ISeqHelper {
  constructor(
    private readonly libHelper: IMonomerLibHelper,
    private readonly helmHelper: IHelmHelper,
    private readonly rdKitModule: RDModule
  ) {}

  getHelmToMolfileConverter(monomerLib: IMonomerLibBase): HelmToMolfileConverter {
    return new HelmToMolfileConverter(this.helmHelper, this.rdKitModule, monomerLib);
  }

  async helmToAtomicLevel(
    helmCol: DG.Column<string>, chiralityEngine?: boolean, highlight?: boolean, overrideMonomerLib?: IMonomerLibBase
  ): Promise<ToAtomicLevelRes> {
    const monomerLib: IMonomerLibBase = overrideMonomerLib ?? this.libHelper.getMonomerLib();

    const df: DG.DataFrame = helmCol.dataFrame;
    const molColName: string = getMolColName(df, helmCol.name);

    const converter = this.getHelmToMolfileConverter(monomerLib);

    //#region From HelmToMolfileConverter.convertToRdKitBeautifiedMolfileColumn

    const molfilesV3K = converter.convertToMolfileV3K(helmCol, this.rdKitModule, monomerLib);

    const beautifiedMolList: (RDMol | null)[] = molfilesV3K.map((item) => {
      const molfile = item.molfile;
      if (molfile === '')
        return null;
      const mol = this.rdKitModule.get_mol(molfile);
      if (!mol)
        return null;
      mol.set_new_coords();
      mol.normalize_depiction(1);
      mol.straighten_depiction(true);
      return mol;
    });

    let molList: string[];
    if (chiralityEngine) {
      molList = converter.getMolV3000ViaOCL(beautifiedMolList, molColName).toList();
      // TODO: Cleanup mol objects
    } else {
      molList = beautifiedMolList.map((mol) => {
        if (mol === null)
          return '';
        const molBlock = mol.get_v3Kmolblock();
        mol!.delete();
        return molBlock;
      });
    }

    //#endregion From HelmToMolfileConverter

    const molHlList = molfilesV3K.map((item: MolfileWithMap) => getMolHighlight(item.monomers.values(), monomerLib));

    const molCol = DG.Column.fromStrings(molColName, molList);
    molCol.semType = DG.SEMTYPE.MOLECULE;
    molCol.meta.units = DG.UNITS.Molecule.MOLBLOCK;
    molCol.setTag(ChemTags.SEQUENCE_SRC_COL, helmCol.name);

    return {molCol: molCol, warnings: []};
  }

  static getInstance(): Promise<SeqHelper> {
    let res = window.$seqHelperPromise;
    if (res == undefined) {
      res = window.$seqHelperPromise = (async () => {
        if (!_package.initialized)
          throw new Error('Bio package is not initialized, call Bio:getSeqHelper');
        const instance = new SeqHelper(
          await MonomerLibManager.getInstance(), await getHelmHelper(), _package.rdKitModule);
        return instance;
      })();
    }
    return res;
  }
}
