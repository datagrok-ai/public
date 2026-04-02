import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {RdKitService} from '../rdkit-service/rdkit-service';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';

export type PharmFamilyId = 'Donor' | 'Acceptor' | 'Hydrophobic' | 'Aromatic' | 'Positive' | 'Negative';
export type PharmFamilySet = {[family in PharmFamilyId]: boolean};
export const PHARM_FAMILY_NAMES: PharmFamilyId[] =
  ['Donor', 'Acceptor', 'Hydrophobic', 'Aromatic', 'Positive', 'Negative'];

const FAMILY_LETTER_TO_ID: {[letter: string]: PharmFamilyId} = {
  'D': 'Donor',
  'A': 'Acceptor',
  'H': 'Hydrophobic',
  'a': 'Aromatic',
  'P': 'Positive',
  'N': 'Negative',
};

export async function runPharmacophoreDetection(
  moleculeCol: DG.Column<string>, familySet: PharmFamilySet,
  featuresDf: DG.DataFrame, rdkitService: RdKitService,
): Promise<DG.DataFrame> {
  const familyCol: DG.Column<string> = featuresDf.getCol('family');
  const smartsCol: DG.Column<string> = featuresDf.getCol('smarts');

  const smarts: {[family: string]: string[]} = {};
  for (let i = 0; i < featuresDf.rowCount; i++) {
    const familyLetter = familyCol.get(i);
    const smartsVal = smartsCol.get(i);
    if (!familyLetter || !smartsVal) continue;
    const familyId = FAMILY_LETTER_TO_ID[familyLetter];
    if (familyId && familySet[familyId]) {
      smarts[familyId] ??= [];
      smarts[familyId].push(smartsVal);
    }
  }

  // Reuse the structural alerts worker — it accepts any {groupName: string[]} map
  const result = await (rdkitService as any).getStructuralAlerts(smarts, moleculeCol.toList());

  const resultDf = DG.DataFrame.create(moleculeCol.length);
  for (const [familyName, boolArray] of result)
    resultDf.columns.addNewBool(resultDf.columns.getUnusedName(familyName)).init((i) => boolArray[i]);

  return resultDf;
}

export async function getPharmacophoreFeaturesByFamilies(
  molecules: DG.Column, familySet: PharmFamilySet,
): Promise<DG.DataFrame | undefined> {
  const rdkitService = await chemCommonRdKit.getRdKitService();
  const featuresDf = await grok.data.loadTable(
    chemCommonRdKit.getRdKitWebRoot() + 'files/pharmacophore-features.csv');
  const progress = DG.TaskBarProgressIndicator.create('Detecting pharmacophore features...');
  try {
    return await runPharmacophoreDetection(molecules, familySet, featuresDf, rdkitService);
  } catch (e) {
    grok.shell.error('Pharmacophore feature detection failed');
    grok.log.error(`Pharmacophore feature detection failed: ${e}`);
  } finally {
    progress.close();
  }
}
