import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full.js';

export async function propertiesWidget(smiles: string) {
  const mol = OCL.Molecule.fromSmiles(smiles);
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
    'Name': ui.wait(async () => {
      const iupacName = await grok.functions.call(`PubChemApi:GetIupacName`, {'smiles': smiles});
      return ui.divText(iupacName ?  iupacName.sval : 'Not found in PubChem');
    }),
  };

  return new DG.Widget(ui.tableFromMap(propertiesMap));
}
  