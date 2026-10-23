import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full';
import {oclMol, renderDescription} from '../utils/chem-common-ocl';
import {_convertMolNotation} from '../utils/convert-notation-utils';
import {PackageFunctions} from '../package';
import {OCLService} from '../open-chem/ocl-service';
import {ChemicalEngine, engineInEffect, isCruxRequested} from '../engine/chem-engine';
import {ChemWasm, cruxModuleSync} from '../engine/crux-module';
import {cruxOrRdKit, cruxOrRdKitAsync} from '../engine/fallback';
import {getCruxService} from '../crux/crux-searches';
import {withinSmilesLimit} from '../utils/chem-common';
import '../../css/chem.css';

/**
 * Crux mode's drug-likeness ("Chemical Engine": Crux; crux-core's decision chem-props): RDKit's QED, a 0..1
 * desirability (Bickerton, Paolini, Besnard, Muresan & Hopkins, Nature Chem. 2012, 4, 90), by crux's `chemQed`
 * (identical to RDKit's values on ChEMBL-10k), with the eight properties it weighs, by `QED.properties`' keys.
 * RDKit mode keeps OpenChemLib's score (an open-ended sum of fragment increments).
 */
export const QED_PROPERTIES: {[key: string]: string} = {
  'MW': 'Molecular weight',
  'ALOGP': 'ALogP (Wildman-Crippen)',
  'HBA': 'H-bond acceptors',
  'HBD': 'H-bond donors',
  'PSA': 'Polar surface area',
  'ROTB': 'Rotatable bonds',
  'AROM': 'Aromatic rings',
  'ALERTS': 'Unwanted groups',
};

/** The column of the whole-table QED in Crux mode. */
export const QED_COLUMN = 'QED';

/**
 * In Crux mode with the Crux module ready, one molecule's QED and its properties ({@link QED_PROPERTIES}) by crux's
 * `chemQed`; null otherwise, and for a molecule Crux cannot read (OpenChemLib's score then shows).
 */
export function cruxQed(molecule: string): {[key: string]: number} | null {
  const m = withinSmilesLimit(molecule);
  if (engineInEffect() !== ChemicalEngine.Crux || m === null)
    return null;
  return cruxOrRdKit('drugLikeness', () => {
    const values = JSON.parse((cruxModuleSync() as unknown as ChemWasm).chemQed([m]));
    return values.QED[0] === null ? null : Object.fromEntries(Object.keys(values).map((k) => [k, values[k][0]]));
  }, () => null);
}

/**
 * The drug-likeness of a column's molecules as a column: in Crux mode QED by crux's `chemQed` in the Crux workers (a
 * molecule Crux cannot read, an empty cell or a SMILES over Chem's length limit gives an empty cell), OpenChemLib's
 * score when Crux fails; in RDKit mode OpenChemLib's score.
 */
export async function getDrugLikenessColumn(molCol: DG.Column<string>): Promise<DG.Column> {
  const ocl = async () => {
    const oclService = new OCLService();
    try {
      const res = await oclService.getDrugLikeness(molCol);
      const name = Object.keys(res)[0];
      return DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, name, res[name]);
    } finally {
      oclService.terminate();
    }
  };
  if (!isCruxRequested())
    return ocl();
  return cruxOrRdKitAsync('drugLikeness', async () => {
    const values = await (await getCruxService()).chemQed(molCol.toList().map(withinSmilesLimit));
    return DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, QED_COLUMN, values[QED_COLUMN]);
  }, ocl);
}

/** The Crux panel: QED, the eight properties it weighs, a button adding the whole table's QED. */
function qedWidget(molString: DG.SemanticValue, qed: {[key: string]: number}): DG.Widget {
  const calcForWholeButton = ui.icons.add(async () => {
    const pi = DG.TaskBarProgressIndicator.create('Calculating QED');
    try {
      const col = await getDrugLikenessColumn(molString.cell.column);
      col.name = molString.cell.dataFrame.columns.getUnusedName(col.name);
      molString.cell.dataFrame.columns.add(col);
    } catch (e) {
      console.error(e);
    } finally {
      pi.close();
    }
  }, 'Calculate QED for the whole table');
  calcForWholeButton.classList.add('chem-drug-likeness-calc-all-button');
  const scoreDiv = ui.divH([ui.label(`QED: ${qed[QED_COLUMN].toFixed(2)}`), calcForWholeButton],
    {classes: 'chem-drug-likeness-widget-score-div'});
  const about = ui.divText('Quantitative estimate of drug-likeness, 0 to 1 (Bickerton et al., Nature Chem. 2012), ' +
    'from these properties:');
  const fmt = (v: number) => Number.isInteger(v) ? `${v}` : v.toFixed(2);
  const table = ui.tableFromMap(Object.fromEntries(Object.entries(QED_PROPERTIES)
    .map(([key, label]) => [label, fmt(qed[key])])));
  const host = ui.divV([scoreDiv, about, table], {classes: 'ui-box'});
  ui.tools.setHoverVisibility(host, [calcForWholeButton]);
  return new DG.Widget(host);
}

const _dlp = new OCL.DruglikenessPredictor();

export function assessDruglikeness(molString: string): [number, OCL.IParameterizedString[]] {
  const mol = oclMol(molString);
  return [_dlp.assessDruglikeness(mol), _dlp.getDetail()];
}

export function drugLikenessWidget(molString: DG.SemanticValue): DG.Widget {
  const rdKitModule = PackageFunctions.getRdKitModule();
  let convertedMolstring: string = '';
  try {
    convertedMolstring = _convertMolNotation(molString.value,
      DG.chem.Notation.Unknown, DG.chem.Notation.Smiles, rdKitModule);
  } catch (e) {
    return new DG.Widget(ui.divText('Molecule is possibly malformed'));
  }
  const qed = cruxQed(convertedMolstring);
  if (qed)
    return qedWidget(molString, qed);
  let score: number;
  let description: OCL.IParameterizedString[];
  try {
    [score, description] = assessDruglikeness(convertedMolstring);
  } catch (e) {
    return new DG.Widget(ui.divText('Could not asses drug likeness'));
  }
  const descriptionHost = renderDescription(description);
  descriptionHost.style.overflow = 'hidden';
  descriptionHost.style.maxHeight = '400px';
  const calcForWholeButton = ui.icons.add(async () => {
    const pi = DG.TaskBarProgressIndicator.create('Calculating drug likeness scores');
    try {
      const oclService = new OCLService();
      const res = await oclService.getDrugLikeness(molString.cell.column);
      oclService.terminate();
      Object.keys(res).forEach((k) => {
        molString.cell.dataFrame.columns.add(DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, k, res[k]));
      });
    } catch (e) {
      console.error(e);
    } finally {
      pi.close();
    }
  }, 'Calculate drug likeness for the whole table');
  const scoreDiv = ui.divH([ui.label(`Score: ${score.toFixed(2)}`), calcForWholeButton],
    {classes: 'chem-drug-likeness-widget-score-div'});
  calcForWholeButton.classList.add('chem-drug-likeness-calc-all-button');

  const host = ui.divV([scoreDiv,
    ui.label(` ${description[0].value}`), descriptionHost], {classes: 'ui-box'});
  ui.tools.setHoverVisibility(host, [calcForWholeButton]);

  return new DG.Widget(host);
}
