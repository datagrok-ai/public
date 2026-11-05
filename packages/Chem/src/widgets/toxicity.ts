import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full';
import {oclMol} from '../utils/chem-common-ocl';
import {getRdKitModule} from '../utils/chem-common-rdkit';
import {_convertMolNotation} from '../utils/convert-notation-utils';
import {OCLService, riskLevels, riskTypes} from '../open-chem/ocl-service';
import {ChemicalEngine, engineInEffect, isCruxRequested, requestedEngine} from '../engine/chem-engine';
import {ChemWasm, cruxModuleSync} from '../engine/crux-module';
import {cruxOrRdKit} from '../engine/fallback';
import {getCruxService} from '../crux/crux-searches';
import {ToxicityAlert} from '../crux/crux-api';
import {drawMoleculeToCanvas} from '../utils/chem-common-rdkit';
import {withinSmilesLimit} from '../utils/chem-common';
import '../../css/chem.css';

/**
 * Crux mode's toxicity ("Chemical Engine": Crux; crux-core's decision chem-props): RDKit's FilterCatalog alerts of the
 * Brenk and NIH catalogs by crux's `chemAlerts` (identical to RDKit's hits on ChEMBL-10k), by its catalog keys, with
 * what each holds. The Structural Alerts panel holds the other sets (PAINS and the ChEMBL rule sets). RDKit mode
 * keeps OpenChemLib's four risks, and the Toxicity Risks menu is OpenChemLib's in both modes.
 */
export const TOXICITY_ALERT_CATALOGS: {[catalog: string]: string} = {
  'Brenk': 'Brenk: 105 unwanted groups, for potential toxicity or unfavourable pharmacokinetics (Brenk et al., ' +
    'ChemMedChem 2008, 3, 435)',
  'NIH': 'NIH: 180 problematic functional groups (Jadhav et al., J. Med. Chem. 2010, 53, 37; Doveston et al., Org. ' +
    'Biomol. Chem. 2015, 13, 859)',
};

/** Why the toxicity alerts are not computed in RDKit mode. */
export const TOXICITY_ALERTS_NEED_CRUX = 'Toxicity alerts use RDKit\'s FilterCatalog, which Chem computes with the ' +
  'Crux engine: set Chem\'s Chemical Engine to Crux';

/** Why the toxicity alerts are not computed in Crux mode after the Crux file failed to load or stopped. */
export const TOXICITY_ALERTS_CRUX_FAILED = 'Toxicity alerts use RDKit\'s FilterCatalog, which Chem computes with the ' +
  'Crux engine, and Crux failed in this session (see the console)';

/** The column of a catalog's alert counts. */
export function alertCountColumnName(catalog: string): string {
  return `${catalog} alerts`;
}

/**
 * In Crux mode with the Crux module ready, one molecule's alerts per catalog ({@link TOXICITY_ALERT_CATALOGS}) by
 * crux's `chemAlerts`, each with its matches' atoms; null otherwise, and for a molecule Crux cannot read.
 */
export function cruxToxicityAlerts(molecule: string): {[catalog: string]: ToxicityAlert[]} | null {
  const m = withinSmilesLimit(molecule);
  if (engineInEffect() !== ChemicalEngine.Crux || m === null)
    return null;
  const catalogs = Object.keys(TOXICITY_ALERT_CATALOGS);
  return cruxOrRdKit('toxicity', () => {
    const res = JSON.parse((cruxModuleSync() as unknown as ChemWasm).chemAlerts([m], catalogs, true));
    return res[catalogs[0]][0] === null ? null : Object.fromEntries(catalogs.map((c) => [c, res[c][0]]));
  }, () => null);
}

/**
 * Each catalog's alert counts for a column's molecules, a column per catalog ({@link alertCountColumnName}): by crux's
 * `chemAlerts` in the Crux workers (a molecule Crux cannot read, an empty cell or a SMILES over Chem's length limit
 * gives an empty cell). Crux mode only: RDKit's FilterCatalog is not in RDKit.js, so RDKit mode throws.
 */
export async function getToxicityAlertsColumns(col: DG.Column<string>, catalogs: string[]): Promise<DG.Column[]> {
  if (!isCruxRequested()) {
    throw new Error(requestedEngine() === ChemicalEngine.Crux ? TOXICITY_ALERTS_CRUX_FAILED :
      TOXICITY_ALERTS_NEED_CRUX);
  }
  const res = await (await getCruxService()).chemAlertNames(col.toList().map(withinSmilesLimit), catalogs);
  return catalogs.map((c) => DG.Column.fromList(DG.COLUMN_TYPE.INT, alertCountColumnName(c),
    res[c].map((hits) => hits === null ? null : hits.length)));
}

/** The bonds of `smiles` (as the Crux module reads it) between two of `atoms`, for a match's highlight. */
function bondsAmong(smiles: string, atoms: number[]): number[] {
  const mol = (cruxModuleSync() as unknown as ChemWasm).get_mol(smiles);
  if (!mol)
    return [];
  try {
    const bonds: {atoms: number[]}[] = JSON.parse(mol.get_json()).molecules[0].bonds ?? [];
    const set = new Set(atoms);
    return bonds.flatMap((b, i) => set.has(b.atoms[0]) && set.has(b.atoms[1]) ? [i] : []);
  } finally {
    mol.delete();
  }
}

/** The Crux panel: per catalog none or its alerts, each drawn with its atoms; a button per catalog adding the whole
 * table's counts. */
function toxicityAlertsWidget(molSemValue: DG.SemanticValue, smiles: string,
  alerts: {[catalog: string]: ToxicityAlert[]}): DG.Widget {
  const host = ui.divV([]);
  for (const [catalog, hits] of Object.entries(alerts)) {
    const calcForWholeButton = ui.icons.add(async () => {
      const pi = DG.TaskBarProgressIndicator.create(`Counting ${catalog} alerts...`);
      try {
        const [col] = await getToxicityAlertsColumns(molSemValue.cell.column, [catalog]);
        col.name = molSemValue.cell.dataFrame.columns.getUnusedName(col.name);
        molSemValue.cell.dataFrame.columns.add(col);
      } catch (e: any) {
        grok.shell.error(e?.message ?? `${e}`);
      } finally {
        pi.close();
      }
    }, `Count ${catalog} alerts for the whole table`);
    calcForWholeButton.classList.add('chem-toxicity-widget-calc-all-button');
    const verdict = ui.divText(hits.length ? `${hits.length} alert${hits.length > 1 ? 's' : ''}` : 'None',
      {classes: 'chem-toxicity-widget-risk-row'});
    verdict.style.color = DG.Color.toHtml(hits.length ? DG.Color.darkRed : DG.Color.darkGreen);
    verdict.style.fontWeight = 'bolder';
    const label = ui.label(catalog);
    ui.tooltip.bind(label, TOXICITY_ALERT_CATALOGS[catalog]);
    const header = ui.divH([label, verdict, calcForWholeButton], {style: {gap: '8px', alignItems: 'center'}});
    ui.tools.setHoverVisibility(header, [calcForWholeButton]);
    const list = ui.div(hits.map((hit) => {
      const canvas = ui.canvas(200, 100);
      const atoms = [...new Set(hit.atoms.flat())];
      try {
        drawMoleculeToCanvas(0, 0, 200, 100, canvas, smiles, null, undefined,
          {atoms, bonds: bondsAmong(smiles, atoms)});
      } catch (e) {
        console.warn(e);
      }
      return ui.div([ui.divText(hit.name), canvas], 'd4-flex-col chem-mol-box');
    }), {classes: 'd4-flex-wrap chem-search-panel-wrapper'});
    host.append(header, list);
  }
  return new DG.Widget(host);
}

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

export async function getToxicityRisksColumns(col: DG.Column,
  toxRisks: {[_ in keyof typeof ToxRiskArgs]?: boolean}): Promise<DG.Column[]> {
  // ids are from 0 to 3 in the order of arguments
  const toxRiskIds: number[] = [];

  Object.keys(toxRisks).forEach((key) => {
    toxRisks[key as keyof typeof ToxRiskArgs] && toxRiskIds.push(ToxRiskArgs[key as keyof typeof ToxRiskArgs]);
  });
  const oclService = new OCLService();
  const risks = await oclService.getChemToxicity(col, toxRiskIds);
  oclService.terminate();
  const cols: DG.Column[] = [];
  const {dataFrame} = col;
  toxRiskIds.forEach((riskId) => {
    const baseName = riskTypes[riskId];
    const riskName = dataFrame ? dataFrame.columns.getUnusedName(baseName) : baseName;
    const toxCol = DG.Column.fromStrings(riskName, risks[riskId].map((risk) => riskLevels[risk]));
    toxCol.meta.colors.setCategorical(riskColorCoding);
    cols.push(toxCol);
  });
  return cols;
}

export function toxicityWidget(molSemValue: DG.SemanticValue): DG.Widget {
  const rdKitModule = getRdKitModule();
  let molStr: string = molSemValue.value;
  try {
    molStr = _convertMolNotation(molStr, DG.chem.Notation.Unknown, DG.chem.Notation.Smiles, rdKitModule);
  } catch (e) {
    return new DG.Widget(ui.divText('Molecule is possibly malformed'));
  }
  const alerts = cruxToxicityAlerts(molStr);
  if (alerts)
    return toxicityAlertsWidget(molSemValue, molStr, alerts);
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
          const col = await getToxicityRisksColumns(molSemValue.cell.column, {[riskName]: true});
          if (col.length)
            molSemValue.cell.dataFrame.columns.add(col[0]);
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
