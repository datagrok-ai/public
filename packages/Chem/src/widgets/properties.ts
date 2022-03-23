import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full.js';
import {oclMol} from '../utils/chem-common-ocl';

export async function propertiesWidget(smiles: string) {
  const mol = oclMol(smiles);
  const formula = mol.getMolecularFormula();
  const molProps = new OCL.MoleculeProperties(mol);

  const propertiesMap = {
    'SMILES': smiles,
    'Formula': formula.formula,
    'MW': formula.absoluteWeight,
    'Number of HBA': molProps.acceptorCount,
    'Number of HBD': molProps.donorCount,
    'LogP': molProps.logP,
    'LogS': molProps.logS,
    'Polar Surface Area': molProps.polarSurfaceArea,
    'Number of rotatabe bonds': molProps.rotatableBondCount,
    'Number of stereo centers': molProps.stereoCenterCount,
    'Name': ui.wait(async () => ui.divText(await getIUPACName(smiles))),
  };

  return new DG.Widget(ui.tableFromMap(propertiesMap));
}

async function getIUPACName(smiles: string): Promise<string> {
  const url = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${smiles}/property/IUPACName/JSON`;
  const response = await fetch(url);
  const responseJson = await response.json();
  const result = responseJson['PropertyTable']['Properties'][0];
  return result.hasOwnProperty('IUPACName') ? result['IUPACName'] : 'Not found in PubChem';
}
