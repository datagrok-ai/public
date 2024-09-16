import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {getAvailableMonomers, getAvailableMonomerMols} from './utils';
import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';

export const PT_CHEM_EXAMPLE = `


 22 24  0  0  0  0  0  0  0  0999 V2000
    0.3128   -0.7509    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3128    0.0740    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4054   -1.1623    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1081   -0.7509    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4054    0.4877    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1081    0.0740    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8175   -1.1623    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0222    0.4877    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8175    0.4877    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5292   -0.7509    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0222    1.3127    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7227    1.7263    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4054   -1.9896    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5292    0.0740    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4544    1.3127    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7406    0.0740    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0222   -1.1623    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.4544    0.4877    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8175   -1.9896    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2453    0.4877    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.1670    1.7285    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7149    2.5513    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  2  1  2  0  0  0  0
  3  1  1  0  0  0  0
  4  3  1  0  0  0  0
  5  2  1  0  0  0  0
  6  5  1  0  0  0  0
  7  4  1  0  0  0  0
  8  2  1  0  0  0  0
  9  6  1  0  0  0  0
 10  7  2  0  0  0  0
 11  8  2  0  0  0  0
 12 11  1  0  0  0  0
 13  3  2  0  0  0  0
 14  9  2  0  0  0  0
 15 18  1  0  0  0  0
 16  8  1  0  0  0  0
 17  1  1  0  0  0  0
 18 16  2  0  0  0  0
 19  7  1  0  0  0  0
 20 14  1  0  0  0  0
  6  4  2  0  0  0  0
 15 12  2  0  0  0  0
 14 10  1  0  0  0  0
 15 21  1  0  0  0  0
 12 22  1  0  0  0  0
M  RGP  1  22   1
M  END`;

export async function getEnumerationChem(molString: string, screenLibrary: string):
  Promise<string[]> {
  const variableMonomers = await getAvailableMonomers(screenLibrary);
  const variableMols = await getAvailableMonomerMols(screenLibrary);
  const enumerations = new Array<string>(variableMonomers.length);

  const rdkitModule: RDModule = await grok.functions.call('Chem:getRdKitModule');
  const molScaffold: RDMol = rdkitModule.get_mol(molString);
  const smiScaffold = molScaffold.get_smiles();
  molScaffold.delete();

  const smilesSubsts = new Array<string>(variableMonomers.length);

  for (let i = 0; i < variableMonomers.length; i++) {
    const name = variableMonomers[i];
    const molBlock = variableMols[name];
    const molSubst: RDMol = rdkitModule.get_mol(molBlock);
    smilesSubsts[i] = molSubst.get_smiles();
    molSubst.delete();
  }

  for (let i = 0; i < variableMonomers.length; i++) {
    let molRes: RDMol | null = null;
    try {
      //TODO: use RDKit linking function when exposed
      const smiResRaw = `${smiScaffold}.${smilesSubsts[i]}`.replaceAll('[1*]C', 'C([1*])').replaceAll('[1*]c', 'c([1*])').replaceAll('[1*]O', 'O([1*])').replaceAll('[1*]N', 'N([1*])');
      const smiRes = `${smiResRaw}`.replaceAll('([1*])', '9').replaceAll('[1*]', '9');
      molRes = rdkitModule.get_mol(smiRes, JSON.stringify({mappedDummiesAreRGroups: true}));
      const molV3 = molRes.get_v3Kmolblock();
      enumerations[i] = molV3;
    } catch (err:any) {
      enumerations[i] = '';
    } finally {
      molRes?.delete();
    }
  }


  return enumerations;
}
