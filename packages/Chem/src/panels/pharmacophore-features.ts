import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {RdKitService} from '../rdkit-service/rdkit-service';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {HIGHLIGHT_BY_SCAFFOLD_TAG} from '../constants';
import {FAMILY_INFO} from '../widgets/pharmacophore-features';
import {IColoredScaffold} from '../rendering/rdkit-cell-renderer';

export type PharmFamilyId = 'Donor' | 'Acceptor' | 'Hydrophobic' | 'Aromatic' | 'Positive' | 'Negative' | 'Halogen Bond';
export type PharmFamilySet = {[family in PharmFamilyId]: boolean};
export const PHARM_FAMILY_NAMES: PharmFamilyId[] =
  ['Donor', 'Acceptor', 'Hydrophobic', 'Aromatic', 'Positive', 'Negative', 'Halogen Bond'];

const FAMILY_LETTER_TO_ID: {[letter: string]: PharmFamilyId} = {
  'D': 'Donor',
  'A': 'Acceptor',
  'H': 'Hydrophobic',
  'a': 'Aromatic',
  'P': 'Positive',
  'N': 'Negative',
  'X': 'Halogen Bond',
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

export function buildPharmacophoreHighlightColumn(
  moleculeCol: DG.Column<string>, familySet: PharmFamilySet, featuresDf: DG.DataFrame,
): DG.Column<string> {
  const familyCol: DG.Column<string> = featuresDf.getCol('family');
  const smartsCol: DG.Column<string> = featuresDf.getCol('smarts');
  const scaffolds: IColoredScaffold[] = [];
  for (let i = 0; i < featuresDf.rowCount; i++) {
    const letter = familyCol.get(i);
    const smartsVal = smartsCol.get(i);
    if (!letter || !smartsVal) continue;
    const familyId = FAMILY_LETTER_TO_ID[letter];
    const info = FAMILY_INFO[letter];
    if (familyId && info && familySet[familyId])
      scaffolds.push({color: info.color, molecule: smartsVal});
  }
  const clone = moleculeCol.clone();
  clone.name = 'Pharmacophore';
  clone.semType = DG.SEMTYPE.MOLECULE;
  clone.setTag(HIGHLIGHT_BY_SCAFFOLD_TAG, JSON.stringify(scaffolds));
  return clone;
}

export async function getPharmacophoreFeaturesByFamilies(
  molecules: DG.Column<string>, familySet: PharmFamilySet,
): Promise<DG.DataFrame | undefined> {
  const rdkitService = await chemCommonRdKit.getRdKitService();
  const featuresDf = await grok.data.loadTable(
    chemCommonRdKit.getRdKitWebRoot() + 'files/pharmacophore-features.csv');
  const progress = DG.TaskBarProgressIndicator.create('Detecting pharmacophore features...');
  try {
    const resultDf = await runPharmacophoreDetection(molecules, familySet, featuresDf, rdkitService);
    resultDf.columns.insert(buildPharmacophoreHighlightColumn(molecules, familySet, featuresDf), 0);
    return resultDf;
  } catch (e) {
    grok.shell.error('Pharmacophore feature detection failed');
    grok.log.error(`Pharmacophore feature detection failed: ${e}`);
  } finally {
    progress.close();
  }
}
