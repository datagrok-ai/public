/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package';
import {ALPHABET, ALIGNMENT, NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {OrgType} from '@datagrok-libraries/bio/src/helm/types';

import {
  IMonomerLibFileManager, IMonomerLibHelper
} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';

declare const org: OrgType;
const LIB_PATH = 'System:AppData/Bio/monomer-libraries/';

export function _setPeptideColumn(col: DG.Column): void {
  addCommonTags(col);
  col.meta.units = NOTATION.SEPARATOR;
  col.setTag('separator', '-');
  // col.setTag('cell.renderer', 'sequence');
}

function addCommonTags(col: DG.Column<any>) {
  col.semType = DG.SEMTYPE.MACROMOLECULE;
  col.setTag('aligned', ALIGNMENT.SEQ);
  col.setTag('alphabet', ALPHABET.PT);
}

export async function getAvailableMonomers(screenLibrary: string): Promise<string[]> {
  const monomerLibHelper: IMonomerLibHelper = await grok.functions.call('Bio:getMonomerLibHelper', {});
  const monomerLib = await monomerLibHelper.readLibrary(LIB_PATH, screenLibrary);
  //NOTICE: works with Peptides only
  return monomerLib.getMonomerSymbolsByType('PEPTIDE');
}

export async function getAvailableMonomerMols(screenLibrary: string): Promise<{[monomerSymbol: string]: string}> {
  const monomerLibHelper: IMonomerLibHelper = await grok.functions.call('Bio:getMonomerLibHelper', {});
  const monomerLib = await monomerLibHelper.readLibrary(LIB_PATH, screenLibrary);
  const monomers = monomerLib.getMonomerSymbolsByType('PEPTIDE');
  const mols = new Array<string>(monomers.length);

  //NOTICE: works with Peptides only
  return monomerLib.getMonomerMolsByPolymerType('PEPTIDE')!;
}

export async function getLibrariesList(): Promise<string[]> {
  const monomerLibHelper: IMonomerLibHelper = await grok.functions.call('Bio:getMonomerLibHelper', {});
  const monomerFileManager: IMonomerLibFileManager = await monomerLibHelper.getFileManager();
  return monomerFileManager.getValidLibraryPaths();
}
