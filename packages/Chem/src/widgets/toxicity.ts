import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full';
import {renderDescription} from '../utils/chem-common-ocl';
import {oclMol} from '../utils/chem-common-ocl';
import {getRdKitModule} from '../utils/chem-common-rdkit';
import {_convertMolNotation} from '../utils/convert-notation-utils';
import {OCLService, riskLevels, riskTypes} from '../OCL-service';

const riskColorCoding: {[index: string]: number} = {
  'Unknown': DG.Color.black,
  'None': DG.Color.darkGreen,
  'Low': DG.Color.orange,
  'High': DG.Color.darkRed,
};

export function getRisks(molStr: string): {[index: string]: string} {
  const mol = oclMol(molStr);
  const toxicityPredictor = new OCL.ToxicityPredictor();

  const risks: {[index: string]: string} = {};
  for (const typeId of Object.keys(riskTypes)) {
    const numTypeId = parseInt(typeId);
    const currentRisk = riskLevels[toxicityPredictor.assessRisk(mol, numTypeId)];
    risks[riskTypes[numTypeId]] = currentRisk;
  }

  return risks;
}

export async function addRisksAsColumns(table: DG.DataFrame, col: DG.Column,
  mutagenicity?: boolean, tumorigenicity?: boolean, irritatingEffects?: boolean, reproductiveEffects?: boolean) {
  // ids are from 0 to 3 in the order of arguments
  const toxRiskIds: number[] = [];
  [mutagenicity, tumorigenicity, irritatingEffects, reproductiveEffects].forEach((val, index) => {
    val && toxRiskIds.push(index);
  });
  const oclService = new OCLService();
  const risks = await oclService.getChemToxicity(col, toxRiskIds);
  oclService.terminate();
  toxRiskIds.forEach((riskId) => {
    const riskName = riskTypes[riskId];
    table.columns.add(DG.Column.fromStrings(riskName, risks[riskId].map((risk) => riskLevels[risk])));
  });
}

export function toxicityWidget(molStr: string): DG.Widget {
  const rdKitModule = getRdKitModule();
  try {
    molStr = _convertMolNotation(molStr, DG.chem.Notation.Unknown, DG.chem.Notation.Smiles, rdKitModule);
  } catch (e) {
    return new DG.Widget(ui.divText('Molecule is possibly malformed'));
  }
  let risks: {[key: string]: string};
  try {
    risks = getRisks(molStr);
  } catch (e) {
    return new DG.Widget(ui.divText('Could not analyze toxicity'));
  }

  const risksTable: {[index: string]: HTMLDivElement} = {};
  for (const [type, risk] of Object.entries(risks)) {
    const currentRiskHost = ui.divText(risk);
    currentRiskHost.style.color = DG.Color.toHtml(riskColorCoding[risk]);
    currentRiskHost.style.fontWeight = 'bolder';

    //@ts-ignore
    ui.tooltip.bind(currentRiskHost, () => renderDescription(toxicityPredictor.getDetail(mol, numTypeId)));

    risksTable[type] = currentRiskHost;
  }

  return new DG.Widget(ui.tableFromMap(risksTable));
}
