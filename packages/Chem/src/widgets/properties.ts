import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full';
import {oclMol} from '../utils/chem-common-ocl';

export function propertiesWidget(smiles: string) {
  const propertiesMap = getPropertiesMap(smiles);
  return new DG.Widget(ui.tableFromMap(propertiesMap));
}

async function getIUPACName(smiles: string): Promise<string> {
  const url = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${smiles}/property/IUPACName/JSON`;
  const response = await fetch(url);
  const responseJson = await response.json();
  const result = responseJson.PropertyTable?.Properties;
  return (result && result[0].hasOwnProperty('IUPACName')) ? result[0].IUPACName : 'Not found in PubChem';
}

export function getPropertiesMap(smiles: string): {
  SMILES: string;
  Formula: string;
  MW: number;
  'Number of HBA': number;
  'Number of HBD': number;
  LogP: number;
  LogS: number;
  'Polar Surface Area': number;
  'Number of rotatable bonds': number;
  'Number of stereo centers': number;
  Name: any;
} {
  const mol = oclMol(smiles);
  const formula = mol.getMolecularFormula();
  const molProps = new OCL.MoleculeProperties(mol);

  return {
    'SMILES': smiles,
    'Formula': formula.formula,
    'MW': formula.absoluteWeight,
    'Number of HBA': molProps.acceptorCount,
    'Number of HBD': molProps.donorCount,
    'LogP': molProps.logP,
    'LogS': molProps.logS,
    'Polar Surface Area': molProps.polarSurfaceArea,
    'Number of rotatable bonds': molProps.rotatableBondCount,
    'Number of stereo centers': molProps.stereoCenterCount,
    'Name': ui.wait(async () => ui.divText(await getIUPACName(smiles))),
  };
}
