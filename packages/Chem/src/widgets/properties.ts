import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full';
import {oclMol} from '../utils/chem-common-ocl';
import {div} from 'datagrok-api/ui';
import $ from 'cash-dom';
import {_convertMolNotation} from '../utils/convert-notation-utils';
import {getRdKitModule} from '../utils/chem-common-rdkit';
import {addCopyIcon} from '../utils/ui-utils';
import {IChemProperty, OCLService, CHEM_PROP_MAP as PROP_MAP} from '../open-chem/ocl-service';
import * as api from '../package-api';

const DGTypeMap = {
  'float': DG.TYPE.FLOAT,
  'int': DG.TYPE.INT,
} as const;

export function getChemPropertyFunc(name: string) : null | ((smiles: string) => any) {
  const p: IChemProperty = PROP_MAP[name];
  if (p !== undefined) {
    return (smiles: string) => {
      const mol = oclMol(smiles);
      return p.valueFunc(mol);
    };
  }
  return null;
}

export function propertiesWidget(semValue: DG.SemanticValue<string>): DG.Widget {
  const rdKitModule = getRdKitModule();
  let smiles: string;
  try {
    smiles = _convertMolNotation(semValue.value, DG.chem.Notation.Unknown,
      DG.chem.Notation.Smiles, rdKitModule);
  } catch (e) {
    return new DG.Widget(ui.divText('Molecule is possibly malformed'));
  }
  const host = div();
  try {
    const propertiesMap = getPropertiesMap(semValue, host);
    const widgetPropMap: {[key: string]: HTMLElement} = {};
    Object.keys(propertiesMap).forEach((it) => widgetPropMap[it] =
      ui.divH([propertiesMap[it].addColumnIcon, propertiesMap[it].value], {style: {'position': 'relative'}}));
    host.appendChild(ui.tableFromMap(widgetPropMap));

    addCopyIcon(widgetPropMap, 'Properties');
  } catch {
    return new DG.Widget(ui.divText('Could not analyze properties'));
  }

  return new DG.Widget(host);
}

export function getPropertiesMap(semValue: DG.SemanticValue<string>, host?: HTMLElement):
 {[k: string]: {value: any, addColumnIcon: HTMLElement | null}} {
  const mol = oclMol(semValue.value);

  function prop(p: IChemProperty, mol: OCL.Molecule): {value: any, addColumnIcon: HTMLElement | null} {
    let addColumnIcon: HTMLElement | null = null;
    if (host && semValue.cell.dataFrame && semValue.cell.column) {
      addColumnIcon = ui.iconFA('plus', async () => {
        const func = DG.Func.find({name: 'addChemPropertiesColumns'})[0];
        if (!func)
          throw new Error('Function addChemPropertiesColumns not found');
        await func.prepare({
          table: semValue.cell.dataFrame,
          molecules: semValue.cell.column,
          MW: p.name === 'MW',
          HBA: p.name === 'HBA',
          HBD: p.name === 'HBD',
          logP: p.name === 'LogP',
          logS: p.name === 'LogS',
          PSA: p.name === 'PSA',
          rotatableBonds: p.name === 'Rotatable bonds',
          stereoCenters: p.name === 'Stereo centers',
          moleculeCharge: p.name === 'Molecule charge',
        }).call(undefined, undefined, {processed: false});
      }, `Calculate ${p.name} for the whole table`);

      ui.tools.setHoverVisibility(host, [addColumnIcon]);
      $(addColumnIcon)
        .css('color', '#2083d5')
        .css('position', 'absolute')
        .css('top', '2px')
        .css('left', '-12px')
        .css('margin-right', '5px');
    }

    return {value: p.valueFunc(mol), addColumnIcon: addColumnIcon};
  }

  const map : {[k: string]: {value: any, addColumnIcon: HTMLElement | null}} = {};
  const props = Object.keys(PROP_MAP);
  for (let n = 0; n < props.length; ++n)
    map[props[n]] = prop(PROP_MAP[props[n]], mol);

  return map;
}

export async function addPropertiesAsColumns(df: DG.DataFrame, smilesCol: DG.Column<string>,
  propsMap: string[]) {
  const resCols = await getPropertiesAsColumns(smilesCol, propsMap);
  resCols.forEach((col) => {
    const origColName = col.name;
    col.name = df.columns.getUnusedName(col.name);
    df.columns.add(col);
    col.setTag('CHEM_WIDGET_PROPERTY', origColName);
    col.setTag('CHEM_ORIG_MOLECULE_COLUMN', smilesCol.name);
  });
  return resCols;
}

export async function getPropertiesAsColumns(smilesCol: DG.Column<string>, propsMap: string[]): Promise<DG.Column[]> {
  const oclService = new OCLService();
  const props = await oclService.getChemProperties(smilesCol, propsMap);
  oclService.terminate();
  return propsMap.map((p) => DG.Column.fromList(DGTypeMap[PROP_MAP[p].type], p, props[p]));
}
