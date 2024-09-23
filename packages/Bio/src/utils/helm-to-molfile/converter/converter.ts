/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full';

import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types/index';
import {IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {getHelmHelper, IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {MolfileWithMap, MonomerMap} from '@datagrok-libraries/bio/src/monomer-works/types';

import {Polymer} from './polymer';
import {GlobalMonomerPositionHandler} from './position-handler';

import {_package} from '../../../package';


export class HelmToMolfileConverter {
  constructor(
    private helmColumn: DG.Column<string>,
    private df: DG.DataFrame,
    private libHelper: IMonomerLibHelper,
    private helmHelper: IHelmHelper,
  ) { }

  convertToSmiles(rdKitModule: RDModule): DG.Column<string> {
    const smiles = this.getSmilesList(rdKitModule);
    const columnName = this.df.columns.getUnusedName(`smiles(${this.helmColumn.name})`);
    return DG.Column.fromStrings(columnName, smiles.map((molecule) => {
      if (molecule === null)
        return '';
      return molecule;
    }));
  }

  private getSmilesList(rdKitModule: RDModule): string[] {
    const molfilesV2K = this.convertToMolfileV3KColumn(rdKitModule).toList();
    const smiles = molfilesV2K.map((mol) => DG.chem.convert(mol, DG.chem.Notation.MolBlock, DG.chem.Notation.Smiles));
    return smiles;
  }

  public getMolV3000ViaOCL(beautifiedMols: (RDMol | null)[], columnName: string): DG.Column<string> {
    const beautifiedMolV2000 = beautifiedMols.map((mol) => {
      if (mol === null)
        return '';
      const molBlock = mol.get_v3Kmolblock();
      mol!.delete();
      return molBlock;
    });
    const molv3000Arr = new Array<string>(beautifiedMolV2000.length);
    const chiralityPb = DG.TaskBarProgressIndicator.create(`Handling chirality...`);
    for (let i = 0; i < beautifiedMolV2000.length; i++) {
      const oclMolecule = OCL.Molecule.fromMolfile(beautifiedMolV2000[i]);
      const molV3000 = oclMolecule.toMolfileV3();
      molv3000Arr[i] = molV3000.replace('STERAC1', 'STEABS');
      const progress = i / beautifiedMolV2000.length * 100;
      chiralityPb.update(progress, `${progress?.toFixed(2)}% of molecules completed`);
    }
    chiralityPb.close();
    return DG.Column.fromStrings(columnName, molv3000Arr);
  }

  // @deprecated Use SeqHelper.helmToAtomicLevel
  convertToRdKitBeautifiedMolfileColumn(chiralityEngine: boolean, rdKitModule: RDModule): DG.Column<string> {
    const molfilesV3K = this.convertToMolfileV3KColumn(rdKitModule).toList();
    const beautifiedMols = molfilesV3K.map((item) => {
      if (item === '')
        return null;
      const mol = rdKitModule.get_mol(item);
      if (!mol)
        return null;
      mol.set_new_coords();
      mol.normalize_depiction(1);
      mol.straighten_depiction(true);
      return mol;
    });
    const columnName = this.df.columns.getUnusedName(`molfile(${this.helmColumn.name})`);

    if (chiralityEngine)
      return this.getMolV3000ViaOCL(beautifiedMols, columnName);
    return DG.Column.fromStrings(columnName, beautifiedMols.map((mol) => {
      if (mol === null)
        return '';
      const molBlock = mol.get_v3Kmolblock();
      mol!.delete();
      return molBlock;
    }));
  }


  public convertToMolfileV3KColumn(rdKitModule: RDModule): DG.Column<string> {
    const polymerGraphColumn: DG.Column<string> = this.getPolymerGraphColumn();
    const monomerLib: IMonomerLib = this.libHelper.getMonomerLib();
    const molfileList = polymerGraphColumn.toList().map(
      (pseudoMolfile: string, idx: number) => {
        const helm = this.helmColumn.get(idx);
        if (!helm) return '';

        let resMolfileWithMap: MolfileWithMap;
        try {
          resMolfileWithMap = this.getPolymerMolfile(helm, pseudoMolfile, rdKitModule, monomerLib);
        } catch (err: any) {
          const [errMsg, errStack] = errInfo(err);
          _package.logger.error(errMsg, undefined, errStack);
          resMolfileWithMap = MolfileWithMap.createEmpty();
        }
        return resMolfileWithMap.molfile;
      });
    const molfileColName = this.df.columns.getUnusedName(`molfileV2K(${this.helmColumn.name})`);
    const molfileColumn = DG.Column.fromList('string', molfileColName, molfileList);
    return molfileColumn;
  }

  /** Gets list of monomer molfiles */
  public convertToMolfileV3K(rdKitModule: RDModule): MolfileWithMap[] {
    const polymerGraphColumn: DG.Column<string> = this.getPolymerGraphColumn();
    const monomerLib: IMonomerLib = this.libHelper.getMonomerLib();
    const rowCount = this.helmColumn.length;
    const resList: MolfileWithMap[] = new Array<MolfileWithMap>(rowCount);
    for (let rowIdx = 0; rowIdx < rowCount; ++rowIdx) {
      const helm = this.helmColumn.get(rowIdx);
      if (!helm) {
        resList[rowIdx] = MolfileWithMap.createEmpty();
        continue;
      }

      const pseudoMolfile = polymerGraphColumn.get(rowIdx)!;
      let resMolfile: MolfileWithMap;
      try {
        resMolfile = this.getPolymerMolfile(helm, pseudoMolfile, rdKitModule, monomerLib);
      } catch (err: any) {
        const [errMsg, errStack] = errInfo(err);
        _package.logger.error(errMsg, undefined, errStack);
        resMolfile = MolfileWithMap.createEmpty();
      }
      resList[rowIdx] = resMolfile;
    }
    return resList;
  }

  private getPolymerGraphColumn(): DG.Column<string> {
    const helmStrList = this.helmColumn.toList();
    const molfileList = this.helmHelper.getMolfiles(helmStrList);
    const molfileCol = DG.Column.fromStrings('mols', molfileList);
    return molfileCol;
  }

  private getPolymerMolfile(
    helm: string,
    polymerGraph: string,
    rdKitModule: RDModule,
    monomerLib: IMonomerLib
  ): MolfileWithMap {
    const woGapsRes = this.helmHelper.removeGaps(helm);
    const woGapsHelm = woGapsRes.resHelm;
    const woGapsReverseMap = new Map<number, number>();
    for (const [orgPosIdx, woGapsPosIdx] of (woGapsRes.monomerMap?.entries() ?? [])) {
      woGapsReverseMap.set(woGapsPosIdx, orgPosIdx);
    }
    const globalPositionHandler = new GlobalMonomerPositionHandler(polymerGraph);
    const woGapsPolymer = new Polymer(woGapsHelm, rdKitModule, monomerLib);
    globalPositionHandler.monomerSymbols.forEach((monomerSymbol: string, monomerIdx: number) => {
      const shift = globalPositionHandler.getMonomerShifts(monomerIdx);
      woGapsPolymer.addMonomer(monomerSymbol, monomerIdx, shift);
    });
    const woGapsMolfile: MolfileWithMap = woGapsPolymer.compileToMolfile();
    const orgMonomerMap = new MonomerMap();
    for (const [woGapsPosIdx, m] of woGapsMolfile.monomers.entries()) {
      const orgPosIdx = woGapsReverseMap.get(woGapsPosIdx)!;
      orgMonomerMap.set(orgPosIdx, m);
    }

    return new MolfileWithMap(woGapsMolfile.molfile, orgMonomerMap);
  }
}
