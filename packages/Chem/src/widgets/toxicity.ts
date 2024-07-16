import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full';
import {oclMol} from '../utils/chem-common-ocl';
import {getRdKitModule} from '../utils/chem-common-rdkit';
import {_convertMolNotation} from '../utils/convert-notation-utils';
import {OCLService, riskLevels, riskTypes} from '../open-chem/ocl-service';
import '../../css/chem.css';

const riskColorCoding: {[index: string]: number} = {
  'Unknown': DG.Color.black,
  'None': DG.Color.darkGreen,
  'Low': DG.Color.orange,
  'High': DG.Color.darkRed,
};

const ToxRiskArgs = {
  mutagenicity: 0,
  tumorigenicity: 1,
  irritatingEffects: 2,
  reproductiveEffects: 3,
} as const;

export function getRisks(molStr: string): {[index: string]: string} {
  const toxicityPredictor = new OCL.ToxicityPredictor();
  const mol = oclMol(molStr);
  const risks: {[index: string]: string} = {};
  for (const typeId of Object.keys(riskTypes)) {
    const numTypeId = parseInt(typeId);
    const currentRisk = riskLevels[toxicityPredictor.assessRisk(mol, numTypeId)];
    risks[riskTypes[numTypeId]] = currentRisk;
  }

  return risks;
}

export async function addRisksAsColumns(table: DG.DataFrame, col: DG.Column,
  toxRisks: {[_ in keyof typeof ToxRiskArgs]?: boolean}) {
  // ids are from 0 to 3 in the order of arguments
  const toxRiskIds: number[] = [];

  Object.keys(toxRisks).forEach((key) => {
    toxRisks[key as keyof typeof ToxRiskArgs] && toxRiskIds.push(ToxRiskArgs[key as keyof typeof ToxRiskArgs]);
  });
  const oclService = new OCLService();
  const risks = await oclService.getChemToxicity(col, toxRiskIds);
  oclService.terminate();
  toxRiskIds.forEach((riskId) => {
    const riskName = riskTypes[riskId];
    const toxCol = DG.Column.fromStrings(riskName, risks[riskId].map((risk) => riskLevels[risk]));
    table.columns.add(toxCol);
    toxCol.meta.colors.setCategorical(riskColorCoding);
  });
}

export function toxicityWidget(molSemValue: DG.SemanticValue): DG.Widget {
  const rdKitModule = getRdKitModule();
  let molStr: string = molSemValue.value;
  try {
    molStr = _convertMolNotation(molStr, DG.chem.Notation.Unknown, DG.chem.Notation.Smiles, rdKitModule);
  } catch (e) {
    return new DG.Widget(ui.divText('Molecule is possibly malformed'));
  }
  try {
    const res = getRisks(molStr);
    const tableDiv = ui.div();
    const risksTable: {[index: string]: HTMLDivElement} = {};
    for (const [type, risk] of Object.entries(res)) {
      const currentRiskHost = ui.divText(risk, {classes: 'chem-toxicity-widget-risk-row'});
      currentRiskHost.style.color = DG.Color.toHtml(riskColorCoding[risk]);
      currentRiskHost.style.fontWeight = 'bolder';
      const numTypeId = parseInt(Object.keys(riskTypes).find((key) => riskTypes[parseInt(key)] === type)!);

      const calcForWholeButton = ui.icons.add(async () => {
        const pi = DG.TaskBarProgressIndicator.create(`Calculating ${type}...`);
        try {
          const riskName = Object.keys(ToxRiskArgs)
            .find((key) => ToxRiskArgs[key as keyof typeof ToxRiskArgs] === numTypeId)!;
          await addRisksAsColumns(molSemValue.cell.dataFrame, molSemValue.cell.column, {[riskName]: true});
        } catch (e) {
          console.error(e);
        } finally {
          pi.close();
        }
      }, `Calculate ${type} for whole table`);
      calcForWholeButton.classList.add('chem-toxicity-widget-calc-all-button');
      currentRiskHost.prepend(calcForWholeButton);
      ui.tools.setHoverVisibility(tableDiv, [calcForWholeButton]);

      risksTable[type] = currentRiskHost;
    }
    tableDiv.append(ui.tableFromMap(risksTable));
    return new DG.Widget(tableDiv);
  } catch (e) {
    return new DG.Widget(ui.divText('Could not analyze toxicity'));
  }
}
