import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import {getRdKitModule, drawMoleculeToCanvas, getRdKitWebRoot} from '../utils/chem-common-rdkit';
import {hexToPercentRgb} from '../utils/chem-common';
import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {ISubstruct} from '@datagrok-libraries/chem-meta/src/types';
import {_convertMolNotation} from '../utils/convert-notation-utils';
import {HIGHLIGHT_BY_SCAFFOLD_TAG} from '../constants';
import {IColoredScaffold} from '../rendering/rdkit-cell-renderer';


let featuresDf: DG.DataFrame | null = null;
const _smartsMap: Map<string, RDMol> = new Map();
let rdKitModule: RDModule | null = null;

const FAMILY_INFO: {[key: string]: {name: string; color: string}} = {
  'D': {name: 'Donor', color: '#2196F3'},
  'A': {name: 'Acceptor', color: '#E53935'},
  'H': {name: 'Hydrophobic', color: '#FFEB3B'},
  'P': {name: 'Positive', color: '#00BCD4'},
  'N': {name: 'Negative', color: '#FF9800'},
  'a': {name: 'Aromatic', color: '#4CAF50'},
  'X': {name: 'Halogen Bond', color: '#9C27B0'},
};

// Priority for overlapping atoms (higher = takes precedence)
const FAMILY_PRIORITY: {[key: string]: number} = {
  'P': 7,
  'N': 6,
  'D': 5,
  'A': 4,
  'X': 3,
  'a': 2,
  'H': 1,
};

export interface PharmFeatureMatch {
  index: number;
  family: string;
  familyName: string;
  featureName: string;
  smarts: string;
  color: string;
  atoms: number[];
  bonds: number[];
}

async function loadPharmDataset(): Promise<void> {
  const path = getRdKitWebRoot() + 'files/pharmacophore-features.csv';
  featuresDf = await grok.data.loadTable(path);
  const smartsCol = featuresDf.getCol('smarts');
  rdKitModule ??= getRdKitModule();

  for (let i = 0; i < smartsCol.length; i++) {
    const currentSmarts = smartsCol.get(i);
    if (!_smartsMap.has(currentSmarts)) {
      try {
        _smartsMap.set(currentSmarts, rdKitModule.get_qmol(currentSmarts));
      } catch (e) {
        console.warn(`Failed to parse pharmacophore SMARTS: ${currentSmarts}`);
      }
    }
  }
}

export async function getPharmacophoreFeatures(molecule: string): Promise<PharmFeatureMatch[]> {
  if (featuresDf == null)
    await loadPharmDataset();
  rdKitModule ??= getRdKitModule();

  const matches: PharmFeatureMatch[] = [];
  let mol: RDMol | null = null;
  try {
    if (!grok.chem.isMolBlock(molecule) && molecule?.length > 5000)
      throw new Error('SMILES string longer than 5000 characters not supported');
    mol = rdKitModule.get_mol(molecule);

    const smartsCol = featuresDf!.getCol('smarts');
    const familyCol = featuresDf!.getCol('family');
    const featureNameCol = featuresDf!.getCol('feature_name');

    for (let i = 0; i < smartsCol.length; i++) {
      const smarts = smartsCol.get(i);
      const subMol = _smartsMap.get(smarts);
      if (!subMol) continue;

      const family = familyCol.get(i);
      const info = family ? FAMILY_INFO[family] : undefined;

      try {
        const matchesJson = mol.get_substruct_matches(subMol);
        if (matchesJson !== '{}') {
          const parsed = JSON.parse(matchesJson);
          // get_substruct_matches returns either a single match object or array of matches
          const matchList = Array.isArray(parsed) ? parsed : [parsed];
          for (const match of matchList) {
            const atoms: number[] = match.atoms ?? [];
            const bonds: number[] = match.bonds ?? [];
            if (atoms.length > 0) {
              matches.push({
                index: i,
                family: family!,
                familyName: info?.name ?? family!,
                featureName: featureNameCol.get(i)!,
                smarts,
                color: info?.color ?? '#999999',
                atoms,
                bonds,
              });
            }
          }
        }
      } catch (e) {
        // Skip invalid matches
      }
    }
  } finally {
    mol?.delete();
  }
  return matches;
}

function buildCombinedSubstruct(matches: PharmFeatureMatch[]): ISubstruct {
  const atomColors: {[key: number]: number[]} = {};
  const bondColors: {[key: number]: number[]} = {};
  const atomPriority: {[key: number]: number} = {};
  const bondPriority: {[key: number]: number} = {};
  const allAtoms: Set<number> = new Set();
  const allBonds: Set<number> = new Set();

  for (const match of matches) {
    const color = hexToPercentRgb(match.color);
    if (!color) continue;
    const priority = FAMILY_PRIORITY[match.family] ?? 0;

    for (const atom of match.atoms) {
      allAtoms.add(atom);
      if (!atomPriority[atom] || priority > atomPriority[atom]) {
        atomColors[atom] = color;
        atomPriority[atom] = priority;
      }
    }
    for (const bond of match.bonds) {
      allBonds.add(bond);
      if (!bondPriority[bond] || priority > bondPriority[bond]) {
        bondColors[bond] = color;
        bondPriority[bond] = priority;
      }
    }
  }

  return {
    atoms: [...allAtoms],
    bonds: [...allBonds],
    highlightAtomColors: atomColors,
    highlightBondColors: bondColors,
  };
}

function createColorLegend(matchedFamilies: Set<string>): HTMLElement {
  const items = Object.entries(FAMILY_INFO)
    .filter(([letter]) => matchedFamilies.has(letter))
    .map(([letter, {name, color}]) => {
      const swatch = ui.div('', {style: {
        width: '12px', height: '12px', borderRadius: '2px',
        backgroundColor: color, display: 'inline-block', marginRight: '6px',
      }});
      const label = ui.div(name, {style: {display: 'inline-block', fontSize: '11px'}});
      return ui.div([swatch, label], {style: {display: 'flex', alignItems: 'center', marginBottom: '2px'}});
    });
  return ui.divV(items, {style: {justifyContent: 'center'}});
}

const NO_HIGHLIGHT = 0;

export async function pharmacophoreFeaturesWidget(molecule: string): Promise<DG.Widget> {
  const colors = [NO_HIGHLIGHT].concat(DG.Color.categoricalPalette.slice(0, 10));
  rdKitModule ??= getRdKitModule();

  let matches: PharmFeatureMatch[] = [];
  try {
    matches = await getPharmacophoreFeatures(molecule);
  } catch (e) {
    console.warn(e);
    return new DG.Widget(ui.divText('Molecule is possibly malformed'));
  }
  if (matches.length === 0)
    return new DG.Widget(ui.divText('No pharmacophore features detected'));

  // Convert SMILES to MolBlock for consistent rendering
  if (!DG.chem.isMolBlock(molecule))
    molecule = _convertMolNotation(molecule, DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock, rdKitModule!);

  // Combined overview image
  const overviewWidth = 255;
  const overviewHeight = 170;
  const overviewCanvas = ui.canvas(overviewWidth, overviewHeight);
  const combinedSubstruct = buildCombinedSubstruct(matches);
  drawMoleculeToCanvas(0, 0, overviewWidth, overviewHeight, overviewCanvas, molecule,
    null, {normalizeDepiction: true, straightenDepiction: true}, combinedSubstruct);

  // Color legend
  const matchedFamilies = new Set(matches.map((m) => m.family));
  const legend = createColorLegend(matchedFamilies);

  // Calculate for whole dataset button
  const calcForWholeButton = ui.button('Calculate for whole dataset', async () => {
    const func = DG.Func.find({package: 'Chem', name: 'pharmacophoreFeaturesTopMenu'})[0];
    const funcCall = func.prepare();
    ui.dialog('Pharmacophore Features')
      .add(await funcCall.getEditor())
      .onOK(() => {
        const args: any = {};
        Object.entries(funcCall.inputs).forEach(([key, val]) => args[key] = val);
        func.apply(args);
      })
      .show({center: true});
  });
  calcForWholeButton.style.justifyContent = 'flex-start';

  // Group matches by family
  const grouped: {[family: string]: PharmFeatureMatch[]} = {};
  for (const match of matches) {
    grouped[match.family] ??= [];
    grouped[match.family].push(match);
  }

  // Build per-family sections
  const thumbWidth = 200;
  const thumbHeight = 100;
  const familyOrder = ['D', 'A', 'H', 'P', 'N', 'a', 'X'];
  const sections: HTMLElement[] = [];

  for (const familyLetter of familyOrder) {
    const familyMatches = grouped[familyLetter];
    if (!familyMatches) continue;

    // Group all matches by feature name — collect all atom/bond indices per feature
    const featureMatchGroups = new Map<string, PharmFeatureMatch[]>();
    for (const m of familyMatches) {
      if (!featureMatchGroups.has(m.featureName))
        featureMatchGroups.set(m.featureName, []);
      featureMatchGroups.get(m.featureName)!.push(m);
    }

    const {name: familyName, color: headerColor} = FAMILY_INFO[familyLetter];
    const matchCount = familyMatches.length;

    const featureItems = [...featureMatchGroups.entries()].map(([featureName, featureMatches]) => {
      const match = featureMatches[0];
      const count = featureMatches.length;
      const description = ui.divText(count > 1 ? `${featureName} (${count})` : featureName);
      const imageHost = ui.canvas(thumbWidth, thumbHeight);

      // Build ISubstruct with all matches highlighted in the family color
      const color = hexToPercentRgb(match.color);
      const allAtoms: number[] = [];
      const allBonds: number[] = [];
      const atomColors: {[key: number]: number[]} = {};
      const bondColors: {[key: number]: number[]} = {};
      for (const fm of featureMatches) {
        for (const a of fm.atoms) {allAtoms.push(a); if (color) atomColors[a] = color;}
        for (const b of fm.bonds) {allBonds.push(b); if (color) bondColors[b] = color;}
      }
      const substruct: ISubstruct = {atoms: allAtoms, bonds: allBonds,
        highlightAtomColors: atomColors, highlightBondColors: bondColors};
      drawMoleculeToCanvas(0, 0, thumbWidth, thumbHeight, imageHost, molecule,
        null, {normalizeDepiction: true, straightenDepiction: true}, substruct);

      const moreBtn = ui.iconFA('ellipsis-v', (e: MouseEvent) => {
        e.stopImmediatePropagation();
        const menu = DG.Menu.popup();
        menu.group('Highlight fragment')
          .items(colors.map((color) => ui.tools.click(getColoredDiv(color), () => {
            const col = grok.shell.tv.dataFrame.currentCol;
            const array: IColoredScaffold[] = col.getTag(HIGHLIGHT_BY_SCAFFOLD_TAG) ?
              JSON.parse(col.getTag(HIGHLIGHT_BY_SCAFFOLD_TAG)) : [];
            const substrIdx = array.findIndex((it) => it.molecule === match.smarts);
            if (substrIdx !== -1) {
              if (color !== NO_HIGHLIGHT)
                array[substrIdx].color = DG.Color.toHtml(color)!;
              else
                array.splice(substrIdx, 1);
            } else {
              if (color !== NO_HIGHLIGHT)
                array.push({molecule: match.smarts, color: DG.Color.toHtml(color)!});
            }
            col.setTag(HIGHLIGHT_BY_SCAFFOLD_TAG, JSON.stringify(array));
            grok.shell.tv.dataFrame?.fireValuesChanged();
          })), () => { });
        menu.show();
      }, 'More');
      $(moreBtn).addClass('chem-mol-view-icon pep-more-icon');

      return ui.div([description,
        ui.divV([moreBtn, imageHost], 'chem-mol-box struct-alerts-mol-box')], 'd4-flex-col');
    });

    const swatch = ui.div('', {style: {
      width: '10px', height: '10px', borderRadius: '2px',
      backgroundColor: headerColor, display: 'inline-block', marginRight: '6px',
    }});
    const header = ui.div([swatch, ui.div(`${familyName} (${matchCount})`, {style: {display: 'inline'}})],
      {style: {fontWeight: 'bold', marginTop: '8px', marginBottom: '4px'}});

    const featureList = ui.div(featureItems, {classes: 'd4-flex-wrap chem-search-panel-wrapper'});
    sections.push(ui.divV([header, featureList]));
  }

  const overviewRow = ui.divH([overviewCanvas, legend],
    {style: {alignItems: 'center', gap: '10px', marginBottom: '8px'}});

  const sectionLabel = ui.divText('Individual pharmacophore features',
    {style: {fontWeight: 'bold', marginBottom: '4px', borderBottom: '1px solid var(--grey-2)'}});

  return new DG.Widget(ui.divV([
    calcForWholeButton,
    overviewRow,
    sectionLabel,
    ...sections,
  ]));
}

function getColoredDiv(color: number): HTMLDivElement {
  return color === NO_HIGHLIGHT ?
    ui.div('None', {style: {width: '100%', minHeight: '20px', marginLeft: '2px'}}) :
    ui.div('', {style:
      {width: '100%', minHeight: '20px', marginRight: '6px', backgroundColor: DG.Color.toHtml(color)}});
}
