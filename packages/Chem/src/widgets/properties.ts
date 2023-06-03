import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full';
import {oclMol} from '../utils/chem-common-ocl';
import {div} from 'datagrok-api/ui';
import $ from 'cash-dom';
import {_convertMolNotation} from '../utils/convert-notation-utils';
import {getRdKitModule} from '../utils/chem-common-rdkit';
import {addCopyIcon} from '../utils/ui-utils';

interface IChemProperty {
  name: string;
  type: DG.ColumnType;
  valueFunc: (mol: OCL.Molecule) => any;
}

const PROP_MAP: any = {
  'MW': {name: 'MW', type: DG.TYPE.FLOAT, valueFunc: (m: OCL.Molecule) => m.getMolecularFormula().absoluteWeight},
  'HBA': {name: 'HBA', type: DG.TYPE.INT, valueFunc: (m: OCL.Molecule) => new OCL.MoleculeProperties(m).acceptorCount},
  'HBD': {name: 'HBD', type: DG.TYPE.INT, valueFunc: (m: OCL.Molecule) => new OCL.MoleculeProperties(m).donorCount},
  'LogP': {name: 'LogP', type: DG.TYPE.FLOAT, valueFunc: (m: OCL.Molecule) => new OCL.MoleculeProperties(m).logP},
  'LogS': {name: 'LogS', type: DG.TYPE.FLOAT, valueFunc: (m: OCL.Molecule) => new OCL.MoleculeProperties(m).logS},
  'PSA': {name: 'PSA', type: DG.TYPE.FLOAT, valueFunc: (m: OCL.Molecule) => new OCL.MoleculeProperties(m).polarSurfaceArea},
  'Rotatable bonds': {name: 'Rotatable bonds', type: DG.TYPE.INT, valueFunc: (m: OCL.Molecule) => new OCL.MoleculeProperties(m).rotatableBondCount},
  'Stereo centers': {name: 'Stereo centers', type: DG.TYPE.INT, valueFunc: (m: OCL.Molecule) => new OCL.MoleculeProperties(m).stereoCenterCount},
  'Molecule charge': {name: 'Molecule charge', type: DG.TYPE.INT, valueFunc: (m: OCL.Molecule) => getMoleculeCharge(m)},
};

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

export function getMoleculeCharge(mol: OCL.Molecule): number {
  const atomsNumber = mol.getAllAtoms();
  let moleculeCharge = 0;
  for (let atomIndx = 0; atomIndx <= atomsNumber; ++atomIndx)
    moleculeCharge += mol.getAtomCharge(atomIndx);

  return moleculeCharge;
}

export function propertiesWidget(semValue: DG.SemanticValue<string>): DG.Widget {
  const rdKitModule = getRdKitModule();
  try {
    semValue.value = _convertMolNotation(semValue.value, DG.chem.Notation.Unknown,
      DG.chem.Notation.Smiles, rdKitModule);
  } catch (e) {
    return new DG.Widget(ui.divText('Molecule is possibly malformed'));
  }
  const host = div();
  let mol: OCL.Molecule | null = null;
  try {
    mol = oclMol(semValue.value);
  } catch {
    return new DG.Widget(ui.divText('Could not analyze properties'));
  }

  function prop(p: IChemProperty, mol: OCL.Molecule) : HTMLElement {
    const addColumnIcon = ui.iconFA('plus', () => {
      const molCol: DG.Column<string> = semValue.cell.column;
      const col : DG.Column = DG.Column.fromType(p.type,
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
    }, `Calculate ${name} for the whole table`);

    ui.tools.setHoverVisibility(host, [addColumnIcon]);
    $(addColumnIcon)
      .css('color', '#2083d5')
      .css('position', 'absolute')
      .css('top', '2px')
      .css('left', '-12px')
      .css('margin-right', '5px');

    return ui.divH([addColumnIcon, p.valueFunc(mol)], {style: {'position': 'relative'}});
  }

  const map = {};
  const props = Object.keys(PROP_MAP);
  for (let n = 0; n < props.length; ++n)
    (map as any)[props[n]] = prop(PROP_MAP[props[n]], mol);

  host.appendChild(ui.tableFromMap(map));

  addCopyIcon(map, 'Properties');

  return new DG.Widget(host);
}
