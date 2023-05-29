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

const CHEM_PROPS : IChemProperty[] = [
  {name: 'MW', type: DG.TYPE.FLOAT, valueFunc: (m: OCL.Molecule) => m.getMolecularFormula().absoluteWeight},
  {name: 'HBA', type: DG.TYPE.INT, valueFunc: (m) => new OCL.MoleculeProperties(m).acceptorCount},
  {name: 'HBD', type: DG.TYPE.INT, valueFunc: (m) => new OCL.MoleculeProperties(m).donorCount},
  {name: 'LogP', type: DG.TYPE.FLOAT, valueFunc: (m) => new OCL.MoleculeProperties(m).logP},
  {name: 'LogS', type: DG.TYPE.FLOAT, valueFunc: (m) => new OCL.MoleculeProperties(m).logS},
  {name: 'PSA', type: DG.TYPE.FLOAT, valueFunc: (m) => new OCL.MoleculeProperties(m).polarSurfaceArea},
  {name: 'Rotatable bonds', type: DG.TYPE.INT, valueFunc: (m) => new OCL.MoleculeProperties(m).rotatableBondCount},
  {name: 'Stereo centers', type: DG.TYPE.INT, valueFunc: (m) => new OCL.MoleculeProperties(m).stereoCenterCount},
  {name: 'Molecule charge', type: DG.TYPE.INT, valueFunc: (m) => getMoleculeCharge(m)},
];

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
    if (CHEM_PROPS[n].name === name) {
      return (smiles: string) => {
        const mol = oclMol(smiles);
        return CHEM_PROPS[n].valueFunc(mol);
      };
    }
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

  function prop(p: IChemProperty, mol: OCL.Molecule) {
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
  for (let n = 0; n < CHEM_PROPS.length; ++n)
    (map as any)[CHEM_PROPS[n].name] = prop(CHEM_PROPS[n], mol);


  host.appendChild(ui.tableFromMap(map));

  addCopyIcon(map, 'Properties');

  return new DG.Widget(host);
}
