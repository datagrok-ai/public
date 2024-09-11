import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {ISeqHelper, ISubstruct, SUBSTRUCT_COL, ToAtomicLevelResType} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {getHelmHelper, IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {getRdKitModule} from '@datagrok-libraries/bio/src/chem/rdkit-module';

import {HelmToMolfileConverter} from '../helm-to-molfile/converter';
import {Column, DataFrame} from 'datagrok-api/dg';
import {MolfileWithMap, MonomerMap} from '../helm-to-molfile/converter/types';

import {_package, getMonomerLibHelper} from '../../package';
import {MonomerLibManager} from '../monomer-lib/lib-manager';

type SeqHelperWindowType = Window & { $seqHelperPromise?: Promise<SeqHelper> };
declare const window: SeqHelperWindowType;

export class SeqHelper implements ISeqHelper {
  constructor(
    private readonly libHelper: IMonomerLibHelper,
    private readonly helmHelper: IHelmHelper,
    private readonly rdKitModule: RDModule
  ) {}

  getHelmToMolfileConverter(df: DataFrame, helmCol: Column<string>) {
    return new HelmToMolfileConverter(helmCol, df, this.libHelper, this.helmHelper);
  }

  async helmToAtomicLevel(
    helmCol: DG.Column<string>, chiralityEngine?: boolean, highlight?: boolean
  ): Promise<ToAtomicLevelResType> {
    const getUnusedName = (df: DG.DataFrame | undefined, colName: string): string => {
      if (!df) return colName;
      return df.columns.getUnusedName(colName);
    };

    const df: DG.DataFrame = helmCol.dataFrame;
    const molColName: string = getUnusedName(df, `molfile(${helmCol})`);
    const molHlColName: string = getUnusedName(df, `~${molColName}-hl`);

    const converter = this.getHelmToMolfileConverter(df, helmCol);

    //#region From HelmToMolfileConverter.convertToRdKitBeautifiedMolfileColumn

    const molfilesV3K = converter.convertToMolfileV3K(this.rdKitModule);

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

    const molHlList = molfilesV3K.map((item: MolfileWithMap) => {
      const mmList: MonomerMap[] = item.monomers;

      const hlAtoms: { [key: number]: number[] } = {};
      const hlBonds: { [key: number]: number[] } = {};

      for (const [mm, mmI] of wu.enumerate(mmList)) {
        if (mmI >= 2) continue;

        const mmColor = [Math.random(), Math.random(), Math.random(), 0.3]; // green color
        for (const mAtom of mm.atoms) {
          hlAtoms[mAtom] = mmColor;
        }
        for (const mBond of mm.bonds) {
          hlBonds[mBond] = mmColor;
        }
      }

      let resSubstruct: ISubstruct = {
        atoms: Object.keys(hlAtoms).map((k) => parseInt(k)),
        bonds: Object.keys(hlBonds).map((k) => parseInt(k)),
        highlightAtomColors: hlAtoms,
        highlightBondColors: hlBonds,
      };

      return resSubstruct;
    });

    const molCol = DG.Column.fromStrings(molColName, molList);
    molCol.semType = 'Molecule';
    molCol.temp[SUBSTRUCT_COL] = molHlColName;
    const molHlCol = DG.Column.fromList(DG.COLUMN_TYPE.OBJECT, molHlColName, molHlList);
    molHlCol.semType = 'Molecule-substruct';

    return {molCol: molCol, molHighlightCol: molHlCol};
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
