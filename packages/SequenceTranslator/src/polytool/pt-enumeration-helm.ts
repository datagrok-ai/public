import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Chain} from './pt-conversion';
import {getAvailableMonomers} from './utils'

export const PT_HELM_EXAMPLE = 'PEPTIDE1{[R].[F].[T].[G].[H].[F].[G].[A].[A].[Y].[P].[E].[NH2]}$$$$';

export async function getEnumerationHelm(helmString: string, helmSelections: number[], screenLibrary: string):
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
