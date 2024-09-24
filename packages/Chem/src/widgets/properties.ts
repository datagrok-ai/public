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
  let mol: OCL.Molecule | null = null;
  try {
    mol = oclMol(smiles);
  } catch {
    return new DG.Widget(ui.divText('Could not analyze properties'));
  }

  function prop(p: IChemProperty, mol: OCL.Molecule) : HTMLElement {
    const addColumnIcon = ui.iconFA('plus', () => {
      const molCol: DG.Column<string> = semValue.cell.column;
      const col : DG.Column = DG.Column.fromType(DGTypeMap[p.type],
        semValue.cell.dataFrame.columns.getUnusedName(p.name), molCol.length)
        .setTag('CHEM_WIDGET_PROPERTY', p.name)
        .setTag('CHEM_ORIG_MOLECULE_COLUMN', molCol.name)
        .init((i) => {
          try {
            if (molCol.isNone(i)) return null;
            const mol = oclMol(molCol.get(i)!);
            return p.valueFunc(mol);
          } catch (_) {
            return null;
          }
        });
      semValue.cell.dataFrame.columns.add(col);
    }, `Calculate ${p.name} for the whole table`);

    ui.tools.setHoverVisibility(host, [addColumnIcon]);
    $(addColumnIcon)
      .css('color', '#2083d5')
      .css('position', 'absolute')
      .css('top', '2px')
      .css('left', '-12px')
      .css('margin-right', '5px');

    return ui.divH([addColumnIcon, p.valueFunc(mol)], {style: {'position': 'relative'}});
  }

  const map : {[k: string]: HTMLElement} = {};
  const props = Object.keys(PROP_MAP);
  for (let n = 0; n < props.length; ++n)
    map[props[n]] = prop(PROP_MAP[props[n]], mol);

  host.appendChild(ui.tableFromMap(map));

  addCopyIcon(map, 'Properties');

  return new DG.Widget(host);
}

export async function addPropertiesAsColumns(df: DG.DataFrame, smilesCol: DG.Column<string>,
  propsMap: string[]) {
  const oclService = new OCLService();
  const props = await oclService.getChemProperties(smilesCol, propsMap);

  oclService.terminate();
  propsMap.forEach((p) => {
    const colName = df.columns.getUnusedName(p);
    const col = DG.Column.fromList(DGTypeMap[PROP_MAP[p].type], colName, props[p]);
    df.columns.add(col);
  });
}

export async function getPropertyForMolecule(molecule: string, prop: string) {  
  const smiles = _convertMolNotation(molecule, DG.chem.Notation.Unknown,
      DG.chem.Notation.Smiles, getRdKitModule());
  const oclMolObj = oclMol(smiles);
  const res = PROP_MAP[prop].valueFunc(oclMolObj);
  return res;
}
