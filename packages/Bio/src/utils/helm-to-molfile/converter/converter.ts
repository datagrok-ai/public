/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as OCL from 'openchemlib/full';

import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {_package} from '../../../package';
import {Polymer} from './polymer';
import {GlobalMonomerPositionHandler} from './position-handler';

export class HelmToMolfileConverter {
  constructor(private helmColumn: DG.Column<string>, private df: DG.DataFrame) {
    this.helmColumn = helmColumn;
  }

  async convertToSmiles(): Promise<DG.Column<string>> {
    const smiles = await this.getSmilesList();
    const columnName = this.df.columns.getUnusedName(`smiles(${this.helmColumn.name})`);
    return DG.Column.fromStrings(columnName, smiles.map((molecule) => {
      if (molecule === null)
        return '';
      return molecule;
    }));
  }

  private async getSmilesList(): Promise<string[]> {
    const molfilesV2K = (await this.convertToMolfileV2KColumn()).toList();
    const smiles = molfilesV2K.map((mol) => DG.chem.convert(mol, DG.chem.Notation.MolBlock, DG.chem.Notation.Smiles));
    return smiles;
  }

  private async getMolV3000ViaOCL(beautifiedMols: (RDMol | null)[], columnName: string) {
    const beautifiedMolV2000 = beautifiedMols.map((mol) => {
      if (mol === null)
        return '';
      const molBlock = mol.get_molblock();
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

  async convertToRdKitBeautifiedMolfileColumn(chiralityEngine?: boolean): Promise<DG.Column<string>> {
    const smiles = await this.getSmilesList();
    const rdKitModule: RDModule = await grok.functions.call('Chem:getRdKitModule');
    const beautifiedMols = smiles.map((item) =>{
      if (item === '')
        return null;
      const mol = rdKitModule.get_mol(item);
      if (!mol)
        return null;
      mol.normalize_depiction(1);
      mol.straighten_depiction(true);
      return mol;
    });
    const columnName = this.df.columns.getUnusedName(`molfile(${this.helmColumn.name})`);

    if (chiralityEngine)
      return await this.getMolV3000ViaOCL(beautifiedMols, columnName);
    return DG.Column.fromStrings(columnName, beautifiedMols.map((mol) => {
      if (mol === null)
        return '';
      const molBlock = mol.get_v3Kmolblock();
      mol!.delete();
      return molBlock;
    }));
  }

  private async convertToMolfileV2KColumn(): Promise<DG.Column<string>> {
    const polymerGraphColumn: DG.Column<string> = await this.getPolymerGraphColumn();
    const rdKitModule = await grok.functions.call('Chem:getRdKitModule');
    const molfileList = polymerGraphColumn.toList().map(
      (pseudoMolfile: string, idx: number) => {
        const helm = this.helmColumn.get(idx);
        if (!helm)
          return '';
        let result = '';
        try {
          result = this.getPolymerMolfile(helm, pseudoMolfile, rdKitModule);
        } catch (err: any) {
          const [errMsg, errStack] = errInfo(err);
          _package.logger.error(errMsg, undefined, errStack);
        } finally {
          return result;
        }
      });
    const molfileColName = this.df.columns.getUnusedName(`molfileV2K(${this.helmColumn.name})`);
    const molfileColumn = DG.Column.fromList('string', molfileColName, molfileList);
    return molfileColumn;
  }

  private async getPolymerGraphColumn(): Promise<DG.Column<string>> {
    const polymerGraphColumn: DG.Column<string> =
      await grok.functions.call('HELM:getMolfiles', {col: this.helmColumn});
    return polymerGraphColumn;
  }

  private getPolymerMolfile(
    helm: string,
    polymerGraph: string,
    rdKitModule: RDModule
  ): string {
    const globalPositionHandler = new GlobalMonomerPositionHandler(polymerGraph);
    const polymer = new Polymer(helm, rdKitModule);
    globalPositionHandler.monomerSymbols.forEach((monomerSymbol: string, monomerIdx: number) => {
      const shift = globalPositionHandler.getMonomerShifts(monomerIdx);
      polymer.addMonomer(monomerSymbol, monomerIdx, shift);
    });
    const polymerMolfile = polymer.compileToMolfile();
    return polymerMolfile;
  }
}

