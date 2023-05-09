import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as OCL from 'openchemlib/full';
import {oclMol} from '../utils/chem-common-ocl';
import {div} from "datagrok-api/ui";
import $ from 'cash-dom';
import {_convertMolNotation} from '../utils/convert-notation-utils';
import {getRdKitModule} from '../utils/chem-common-rdkit';
import { MOL_FORMAT } from '../constants';
import { addCopyIcon } from '../utils/ui-utils';
import { StringUtils } from '@datagrok-libraries/utils/src/string-utils';
import { renderMolecule } from '../rendering/render-molecule';
import { renderMultipleHistograms } from '../../../../js-api/dg';
import { histogram } from '../../../../js-api/src/api/ddt.api.g';
import { findMCS } from '../scripts-api';

interface IChemProperty {
  name: string;
  type: DG.ColumnType;
  valueFunc: (mol: OCL.Molecule) => any;
}

const CHEM_PROPS : IChemProperty[] = [
  {name: 'MW', type: DG.TYPE.FLOAT, valueFunc: (m: OCL.Molecule) => m.getMolecularFormula().absoluteWeight},
  {name: 'HBA', type: DG.TYPE.INT, valueFunc: (m) => new OCL.MoleculeProperties(m).acceptorCount},
  {name: 'HBD', type: DG.TYPE.INT, valueFunc: (m) => new OCL.MoleculeProperties(m).donorCount},
  {name: 'LogP', type: DG.TYPE.FLOAT, valueFunc: (m) => new OCL.MoleculeProperties(m).logP},
  {name: 'LogS', type: DG.TYPE.FLOAT, valueFunc: (m) => new OCL.MoleculeProperties(m).logS},
  {name: 'PSA', type: DG.TYPE.FLOAT, valueFunc: (m) => new OCL.MoleculeProperties(m).polarSurfaceArea},
  {name: 'RotB', type: DG.TYPE.INT, valueFunc: (m) => new OCL.MoleculeProperties(m).rotatableBondCount},
  {name: 'StereoC', type: DG.TYPE.INT, valueFunc: (m) => new OCL.MoleculeProperties(m).stereoCenterCount},
  {name: 'Charge', type: DG.TYPE.INT, valueFunc: (m) => getMoleculeCharge(m)}
]

const MAX_COL_LENGTH = 1000;

export function calcChemProperty(name: string, smiles: string) : any {
  const mol = oclMol(smiles);
  for (let n = 0; n < CHEM_PROPS.length; ++n) {
    if (CHEM_PROPS[n].name === name)
      return CHEM_PROPS[n].valueFunc(mol);
  }
  return null;
}

export function getChemPropertyFunc(name: string) : null | ((smiles: string) => any) {
  for (let n = 0; n < CHEM_PROPS.length; ++n) {
    if (CHEM_PROPS[n].name === name)
      return (smiles: string) => {
      const mol = oclMol(smiles);
      return CHEM_PROPS[n].valueFunc(mol);
    };
  }
  return null;
}

function getPropertyValue(molCol: DG.Column, idx: number, p: IChemProperty): any {
  try {
    if (molCol.isNone(idx)) return null;
    const mol = oclMol(molCol.get(idx)!);
    return p.valueFunc(mol);
  }
  catch (_) {
    return null;
  }
}

function getIcon(p: IChemProperty, molCol: DG.Column): HTMLElement {
  const addColumnIcon = ui.iconFA('plus', () => {
    const col : DG.Column = DG.Column.fromType(p.type, molCol.dataFrame.columns.getUnusedName(p.name), molCol.length)
    .setTag('CHEM_WIDGET_PROPERTY', p.name)
    .setTag('CHEM_ORIG_MOLECULE_COLUMN', molCol.name)
    .init((i) => {
      return getPropertyValue(molCol, i, p);
    });
    molCol.dataFrame.columns.add(col);
  }, `Calculate ${p.name} for the whole table`);
  return addColumnIcon;
}

export function getMoleculeCharge(mol: OCL.Molecule): number {
  const atomsNumber = mol.getAllAtoms();
  let moleculeCharge = 0;
  for (let atomIndx = 0; atomIndx <= atomsNumber; ++atomIndx) {
    moleculeCharge += mol.getAtomCharge(atomIndx);
  }
  return moleculeCharge;
}

export async function statsWidget(molCol: DG.Column<string>): Promise<DG.Widget> {
  const host = ui.div();
  host.style.marginLeft = '-25px';
  const size = Math.min(molCol.length, MAX_COL_LENGTH);
  const randomIndexes = Array.from({length: size}, (_, i) => i);

  function getAverage(p: IChemProperty): HTMLDivElement {
    const propertiesArray = new Array(size);
    for (let idx = 0; idx < randomIndexes.length; ++idx) 
      propertiesArray[idx] = getPropertyValue(molCol, randomIndexes[idx], p);
    const col = DG.Column.fromList(p.type, molCol.dataFrame.columns.getUnusedName(p.name), propertiesArray);

    const addColumnIcon = getIcon(p, molCol); 

    ui.tools.setHoverVisibility(host, [addColumnIcon]);
    $(addColumnIcon).addClass('chem-plus-icon');
    addColumnIcon.style.top = '10px';

    const bins = DG.BitSet.create(col.length);
    bins.setAll(true);
    const hist = histogram(col, bins, true, {bins: 20, logScale: false});
  
    const histRoot = ui.div();
    histRoot.classList.add('chem-distribution-hist');
  
    const hw = 100;
    const hh = 35;
    const m = 4;

    const canvas = ui.canvas(hw, hh);
    const ctx = canvas.getContext('2d', {willReadFrequently : true});
    
    if (ctx) 
      renderMultipleHistograms(ctx, new DG.Rect(m, m, hw - 2 * m, hh - 2 * m), [hist], {fill: false, markerSize: 1});
    histRoot.appendChild(canvas);
    
    const colAvg = StringUtils.formatNumber(col.stats.avg);
    const colStdev = StringUtils.formatNumber(col.stats.stdev);
    const textStats = ui.divText(`${colAvg} ± ${colStdev}`);
    textStats.classList.add('chem-stats-text');
    return ui.divH([
      addColumnIcon, 
      textStats,
      histRoot
    ], { style: {'position': 'relative'}});
  
  };

  let map: { [_: string]: HTMLElement }  = {};
  for (let n = 0; n < CHEM_PROPS.length; ++n) {
    map[CHEM_PROPS[n].name] = getAverage(CHEM_PROPS[n]);
  }

  if (molCol.length === size) {
    await findMCS(molCol.name, molCol.dataFrame, true, true).then((res: string) => {
      map['MCS'] = renderMolecule(res, {renderer: 'RDKit'})
    });
  } else
    host.appendChild(ui.divText('Main statistics for sampled data'));

  host.appendChild(ui.tableFromMap(map));

  return new DG.Widget(host);
}

export function propertiesWidget(semValue: DG.SemanticValue<string>): DG.Widget {
  const rdKitModule = getRdKitModule();
  try {
    semValue.value = _convertMolNotation(semValue.value, 'unknown', MOL_FORMAT.SMILES, rdKitModule);
  } catch (e) {
    return new DG.Widget(ui.divText('Molecule is possibly malformed'));
  }
  let host = div();
  try {
    var mol = oclMol(semValue.value);
  } catch {
    return new DG.Widget(ui.divText('Could not analyze properties'));
  }

  function prop(p: IChemProperty) {
    const molCol: DG.Column<string> = semValue.cell.column;
    const addColumnIcon = getIcon(p, molCol);
    ui.tools.setHoverVisibility(host, [addColumnIcon]);
    addColumnIcon.classList.add('chem-plus-icon');

    return ui.divH([addColumnIcon, p.valueFunc(mol)], { style: {'position': 'relative'}});
  }

  let map  = {};
  for (let n = 0; n < CHEM_PROPS.length; ++n) {
    (map as any)[CHEM_PROPS[n].name] = prop(CHEM_PROPS[n]);
  }

  host.appendChild(ui.tableFromMap(map));

  addCopyIcon(map, 'Properties');

  return new DG.Widget(host);
}