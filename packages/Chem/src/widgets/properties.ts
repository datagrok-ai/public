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
import {prop} from '@datagrok-libraries/utils/src/add-icon-utils';

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

export function getPropertiesMap(
  semValue: DG.SemanticValue<string>,
  host?: HTMLElement,
): { [k: string]: { value: any, addColumnIcon: HTMLElement | null } } {
  const mol = oclMol(semValue.value);
  const map: { [k: string]: { value: any, addColumnIcon: HTMLElement | null } } = {};

  for (const key of Object.keys(PROP_MAP)) {
    const chemProp = PROP_MAP[key];
    const {value, addColumnIcon} = prop(
      () => chemProp.valueFunc(mol),
      host && semValue?.cell?.dart &&
        semValue.cell.dataFrame && semValue.cell.column ? host : null,
      `Calculate ${chemProp.name} for the whole table`,
      async () => {
        if (!semValue?.cell?.dataFrame || !semValue.cell.column) return;
        const res = await addPropertiesAsColumns(
          semValue.cell.dataFrame,
          semValue.cell.column,
          [chemProp.name],
        );
        const col = res[0];
        col.setTag('CHEM_WIDGET_PROPERTY', chemProp.name);
        col.setTag('CHEM_ORIG_MOLECULE_COLUMN', semValue.cell.column.name);
      },
    );
    map[key] = {value, addColumnIcon};
  }

  return map;
}

export async function addPropertiesAsColumns(df: DG.DataFrame, smilesCol: DG.Column<string>,
  propsMap: string[]) {
  const resCols = await getPropertiesAsColumns(smilesCol, propsMap);
  resCols.forEach((col) => {
    col.name = df.columns.getUnusedName(col.name);
    df.columns.add(col);
  });
  return resCols;
}

export async function getPropertiesAsColumns(smilesCol: DG.Column<string>, propsMap: string[]): Promise<DG.Column[]> {
  const oclService = new OCLService();
  const props = await oclService.getChemProperties(smilesCol, propsMap);
  oclService.terminate();
  return propsMap.map((p) => DG.Column.fromList(DGTypeMap[PROP_MAP[p].type], p, props[p]));
}
