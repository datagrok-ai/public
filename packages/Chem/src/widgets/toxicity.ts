import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full.js';
import {renderDescription} from '../chem_common_ocl';
import {oclMol} from './drug-likeness';

export function toxicityWidget(smiles: string) {
  const mol = oclMol(smiles);
  const riskTypes: {[index: number]: string} = {
    0: 'Mutagenicity',
    1: 'Tumorigenicity',
    2: 'Irritating effects',
    3: 'Reproductive effects',
  };
  const riskLevels: {[index: number]: string} = {
    0: 'Unknown',
    1: 'None',
    2: 'Low',
    3: 'High',
  };

  const riskColorCoding: {[index: string]: number} = {
    'Unknown': DG.Color.black,
    'None': DG.Color.darkGreen,
    'Low': DG.Color.orange,
    'High': DG.Color.darkRed,
  };

  const risks: {[index: string]: HTMLDivElement} = {};
  const toxicityPredictor = new OCL.ToxicityPredictor();
  Object.keys(riskTypes).forEach((typeId) => {
    //@ts-ignore
    const currentRisk = riskLevels[toxicityPredictor.assessRisk(mol, typeId)];
    const currentRiskHost = ui.divText(currentRisk);

    currentRiskHost.style.color = DG.Color.toHtml(riskColorCoding[currentRisk]);
    currentRiskHost.style.fontWeight = 'bolder';
    //@ts-ignore
    ui.tooltip.bind(currentRiskHost, () => renderDescription(toxicityPredictor.getDetail(mol, typeId)));

    //@ts-ignore
    risks[riskTypes[typeId]] = currentRiskHost;
  });

  return new DG.Widget(ui.tableFromMap(risks));
}
