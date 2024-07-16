import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {OrgType} from '@datagrok-libraries/bio/src/helm/types';
import {
  IMonomerLibFileManager, IMonomerLibHelper
} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';

import {Chain} from './pt-conversion';

declare const org: OrgType;

const LIB_PATH = 'System:AppData/Bio/monomer-libraries/';
export const PT_HELM_EXAMPLE = 'PEPTIDE1{[R].[F].[T].[G].[H].[F].[G].[A].[A].[Y].[P].[E].[NH2]}$$$$';

export async function getLibrariesList(): Promise<string[]> {
  const monomerLibHelper: IMonomerLibHelper = await grok.functions.call('Bio:getMonomerLibHelper', {});
  const monomerFileManager: IMonomerLibFileManager = await monomerLibHelper.getFileManager();
  return monomerFileManager.getValidLibraryPaths();
}

export async function getEnumeration(helmString: string, helmSelections: number[], screenLibrary: string):
  Promise<string[]> {
  const variableMonomers = await getAvailableMonomers(screenLibrary);
  const chain: Chain = Chain.fromHelm(helmString);
  const size = helmSelections.length * variableMonomers.length;
  const enumerations = new Array<string>(size);

  for (let i = 0; i < helmSelections.length; i++) {
    for (let j = 0; j < variableMonomers.length; j++)
      enumerations[i * variableMonomers.length + j] = chain.getHelmChanged(helmSelections[i], variableMonomers[j]);
  }

  return enumerations;
}

async function getAvailableMonomers(screenLibrary: string): Promise<string[]> {
  const monomerLibHelper: IMonomerLibHelper = await grok.functions.call('Bio:getMonomerLibHelper', {});
  const monomerLib = await monomerLibHelper.readLibrary(LIB_PATH, screenLibrary);
  //NOTICE: works with Peptides only
  return monomerLib.getMonomerSymbolsByType('PEPTIDE');
}
