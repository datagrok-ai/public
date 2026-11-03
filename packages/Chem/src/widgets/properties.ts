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
import {ChemicalEngine, byEngineAsync, engineInEffect} from '../engine/chem-engine';
import {ChemWasm, cruxModuleSync} from '../engine/crux-module';
import {cruxOrRdKit} from '../engine/fallback';
import {getCruxService} from '../crux/crux-searches';
import {withinCruxReader} from '../crux/crux-service';

const DGTypeMap = {
  'float': DG.TYPE.FLOAT,
  'int': DG.TYPE.INT,
  'string': DG.TYPE.STRING,
} as const;

/**
 * How Crux mode ("Chemical Engine": Crux) computes each property: crux's `chemProperties`, RDKit's definitions
 * (identical to RDKit's values on ChEMBL-10k; crux-core's decision chem-props), each with its method and reference.
 * RDKit mode computes them with OpenChemLib, as before (Crux mode's drug-likeness and toxicity:
 * widgets/drug-likeness.ts, widgets/toxicity.ts).
 */
export const CRUX_PROPERTY_METHODS: {[name: string]: string} = {
  'MW': 'Average molecular weight from the IUPAC standard atomic weights (RDKit MolWt). OpenChemLib\'s MW ' +
    '(RDKit mode) is the monoisotopic mass. Reference: IUPAC standard atomic weights, as RDKit tabulates them',
  'HBA': 'Lipinski\'s acceptors: N and O atoms (RDKit NumLipinskiHBA). Reference: Lipinski, Lombardo, Dominy & ' +
    'Feeney, Adv. Drug Deliv. Rev. 1997, 23, 3',
  'HBD': 'Donor atoms: N, O and S bearing an H (RDKit NumHDonors). Reference: RDKit\'s donor SMARTS (Lipinski.cpp)',
  'LogP': 'Wildman-Crippen atom contributions (RDKit MolLogP). Reference: Wildman & Crippen, J. Chem. Inf. ' +
    'Comput. Sci. 1999, 39, 868',
  'LogS': 'ESOL: 0.26121 - 0.74167 LogP - 0.0066139 MW + 0.0034515 rotatable bonds - 0.42625 aromatic proportion, ' +
    'on RDKit\'s descriptors. Reference: Delaney, J. Chem. Inf. Comput. Sci. 2004, 44, 1000; the coefficients ' +
    'refit on RDKit by P. Walters, Practical Cheminformatics, 2018',
  'PSA': 'Topological polar surface area with S and P (RDKit TPSA, includeSandP). Reference: Ertl, Rohde & ' +
    'Selzer, J. Med. Chem. 2000, 43, 3714',
  'Rotatable bonds': 'RDKit\'s strict definition: non-ring single bonds between non-terminal atoms, but not an amide ' +
    'C-N, a bond to a CF3, CCl3, CBr3 or tert-butyl carbon, or an atom of a triple bond (RDKit NumRotatableBonds). ' +
    'Reference: RDKit Lipinski.cpp, NumRotatableBondsOptions::Strict',
  'Stereo centers': 'Tetrahedral stereocentres, assigned or not (RDKit CalcNumAtomStereoCenters). Reference: RDKit ' +
    'assignStereochemistry, flagPossibleStereoCenters',
  'Molecule charge': 'The sum of the atoms\' formal charges (RDKit GetFormalCharge)',
  'Molecular formula': 'Hill order: C, H, then the other elements alphabetically; the charge last (RDKit ' +
    'CalcMolFormula). Reference: Hill, J. Am. Chem. Soc. 1900, 22, 478',
};

/**
 * In Crux mode with the Crux module ready, one molecule's properties by crux's `chemProperties`; null otherwise, and
 * for a molecule Crux cannot read (the caller's OpenChemLib path then gives them).
 */
function cruxProperties(molecule: string, names: string[]): {[name: string]: any} | null {
  const m = withinCruxReader(molecule);
  if (engineInEffect() !== ChemicalEngine.Crux || m === null)
    return null;
  return cruxOrRdKit('properties', () => {
    const values = JSON.parse((cruxModuleSync() as unknown as ChemWasm).chemProperties([m], names));
    return values[names[0]][0] === null ? null : Object.fromEntries(names.map((name) => [name, values[name][0]]));
  }, () => null);
}

export function getChemPropertyFunc(name: string) : null | ((smiles: string) => any) {
  const p: IChemProperty = PROP_MAP[name];
  if (p !== undefined) {
    return (smiles: string) => {
      const crux = cruxProperties(smiles, [name]);
      if (crux)
        return crux[name];
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
    const table = ui.tableFromMap(widgetPropMap);
    if (Object.values(propertiesMap).some((p) => p.method)) { // Crux's: each name tells its method and reference
      for (const row of Array.from(table.querySelectorAll('tr'))) {
        const method = propertiesMap[row.cells[0]?.textContent ?? '']?.method;
        if (method)
          ui.tooltip.bind(row.cells[0], method);
      }
    }
    host.appendChild(table);

    addCopyIcon(widgetPropMap, 'Properties');
  } catch {
    return new DG.Widget(ui.divText('Could not analyze properties'));
  }

  return new DG.Widget(host);
}

/**
 * The properties of one molecule: in Crux mode crux's (`method`: how and after whom, {@link CRUX_PROPERTY_METHODS}),
 * else, and for a molecule Crux cannot read, OpenChemLib's.
 */
export function getPropertiesMap(semValue: DG.SemanticValue<string>, host?: HTMLElement):
 {[k: string]: {value: any, addColumnIcon: HTMLElement | null, method?: string}} {
  const props = Object.keys(PROP_MAP);
  const crux = cruxProperties(semValue.value, props);
  const mol = crux ? null : oclMol(semValue.value);

  function prop(p: IChemProperty, mol: OCL.Molecule | null):
    {value: any, addColumnIcon: HTMLElement | null, method?: string} {
    let addColumnIcon: HTMLElement | null = null;
    if (host && semValue?.cell?.dart && semValue.cell.dataFrame && semValue.cell.column) {
      addColumnIcon = ui.iconFA('plus', async () => {
        const res = await addPropertiesAsColumns(semValue.cell.dataFrame, semValue.cell.column, [p.name]);
        const col = res[0];
        col.setTag('CHEM_WIDGET_PROPERTY', p.name);
        col.setTag('CHEM_ORIG_MOLECULE_COLUMN', semValue.cell.column.name);
      }, `Calculate ${p.name} for the whole table`);

      ui.tools.setHoverVisibility(host, [addColumnIcon]);
      $(addColumnIcon)
        .css('color', '#2083d5')
        .css('position', 'absolute')
        .css('top', '2px')
        .css('left', '-12px')
        .css('margin-right', '5px');
    }

    if (crux)
      return {value: crux[p.name], addColumnIcon: addColumnIcon, method: CRUX_PROPERTY_METHODS[p.name]};
    return {value: p.valueFunc(mol!), addColumnIcon: addColumnIcon};
  }

  const map : {[k: string]: {value: any, addColumnIcon: HTMLElement | null, method?: string}} = {};
  for (let n = 0; n < props.length; ++n)
    map[props[n]] = prop(PROP_MAP[props[n]], mol);

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

/**
 * The properties of a column's molecules as columns (the names, types and order asked for): in Crux mode by crux's
 * `chemProperties` in the Crux workers ({@link CRUX_PROPERTY_METHODS}; a molecule Crux cannot read, an empty cell or a
 * SMILES over Chem's length limit gives empty cells, so that a column holds one method's values), OpenChemLib's
 * workers when Crux fails; in RDKit mode OpenChemLib's.
 */
export async function getPropertiesAsColumns(smilesCol: DG.Column<string>, propsMap: string[]): Promise<DG.Column[]> {
  const ocl = async () => {
    const oclService = new OCLService();
    try {
      return await oclService.getChemProperties(smilesCol, propsMap);
    } finally {
      oclService.terminate();
    }
  };
  const props = await byEngineAsync('properties', async () => {
    return (await getCruxService()).chemProperties(smilesCol.toList().map(withinCruxReader), propsMap);
  }, ocl);
  return propsMap.map((p) => DG.Column.fromList(DGTypeMap[PROP_MAP[p].type], p, props[p]));
}
