import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {closeAllAccordionPanes, demoScaffold, getAccordionPane, openMoleculeDataset,
  openSketcher, scrollTable} from '../utils/demo-utils';
import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';
import {awaitCheck} from '@datagrok-libraries/utils/src/test';
import {_package, PackageFunctions} from '../package';
import {rGroupAnalysis} from '../analysis/r-group-analysis';
import {CLIFFS_DF_NAME, activityCliffsIdx} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import {BitArrayMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/types';
import {ScaffoldTreeViewer} from '../widgets/scaffold-tree';
import {MatchedMolecularPairsViewer} from '../analysis/molecular-matched-pairs/mmp-viewer/mmp-viewer';
import {dockSarMatrixTabs} from '../analysis/sar-matrix/sar-matrix-viewer';
import {CELL_H, CELL_W, COL_HEADER_H, CORE_W} from '../analysis/sar-matrix/sar-matrix-ui-common';


export async function _demoChemOverview(): Promise<void> {
  const sketcherType = DG.chem.currentSketcherType;
  DG.chem.currentSketcherType = 'OpenChemLib';

  const firstCols = [
    'smiles',
    'MolWt',
    'ExactMolWt',
    'NOCount',
    'RingCount',
  ];
  const lastCols = [
    'NumRadicalElectrons',
    'MinPartialCharge',
    'MaxAbsPartialCharge',
    'NHOHCount',
    'NumSaturatedCarbocycles',
    'NumAliphaticHeterocycles',
    'FpDensityMorgan1',
    'NumAromaticHeterocycles',
    'NumValenceElectrons',
    'NumRotatableBonds',
    'NumAromaticCarbocycles',
    'NumAliphaticCarbocycles',
    'NumHDonors',
    'FpDensityMorgan3',
    'NumAromaticRings',
    'HeavyAtomMolWt',
    'NumSaturatedRings',
    'NumHAcceptors',
    'NumHeteroatoms',
    'NumSaturatedHeterocycles',
    'NumAliphaticRings',
    'MaxPartialCharge',
    'FpDensityMorgan2',
    'FractionCSP3',
    'HeavyAtomCount'];

  const demoScript = new DemoScript('Overview', 'Overview of Cheminformatics functionality',
    undefined, {autoStartFirstStep: true, path: 'Cheminformatics/Overview'});
  let table: DG.DataFrame;
  let tv: DG.TableView;
  let propPanel: Element;
  let canvas: HTMLCanvasElement;
  let filters: DG.FilterGroup;
  await demoScript
    .step('Load molecules', async () => {
      tv = await openMoleculeDataset('demo_files/demo_smiles.csv');
      const layoutString = await _package.files.readAsText('demo_files/Overview_demo.layout');
      const layout = DG.ViewLayout.fromJson(layoutString);
      await DG.delay(100);
      tv.loadLayout(layout);
      grok.shell.windows.showHelp = false;
      grok.shell.windows.context.visible = true;
      table = tv.dataFrame;
    }, {description: 'Load dataset with molecule columns', delay: 3000})
    .step('Calculate molecule properties', async () => {
      if (!grok.shell.windows.context.visible)
        grok.shell.windows.context.visible = true;
      await awaitCheck(() => document.getElementsByClassName('grok-entity-prop-panel').length > 0,
        'Property panel did not appear', 5000);
      propPanel = document.getElementsByClassName('grok-entity-prop-panel')[0];
      closeAllAccordionPanes(propPanel!);
      const molColumnName = table.columns.bySemType(DG.SEMTYPE.MOLECULE)!.name;
      table.currentCell = table.cell(0, molColumnName);
    }, {description: `Open any pane on the context panel on the right to calculate corresponding properties. For instance, open 'Structure' -> '3D Structure' and 'Biology' -> 'Toxicity'.
      Click any other molecule in dataset to re-calculate properties.
      `, delay: 3000})
    .step('Fast rendering', async () => {
      await DG.delay(1000);
      canvas = tv.grid.root.getElementsByTagName('canvas')[2];
      await scrollTable(canvas, 20000, 50, 20);
    }, {description: 'Molecules are rendered immediately when scrolling dataset', delay: 2000})
    .step('Filter molecules by substructure', async () => {
      await DG.delay(1000);
      filters = tv.getFiltersGroup();
      await DG.delay(1000);
      const sketcherDlg = await openSketcher(filters.root, 'sketch-link');
      const sketcherInput = sketcherDlg!
        .getElementsByClassName('grok-sketcher-input')[0]?.children[0] as HTMLInputElement;
      sketcherInput.value = 'C1CCCCC1';
      sketcherInput.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter'}));
      await DG.delay(3000);
      Array.from(sketcherDlg!.getElementsByTagName('span')).find((el) => el.textContent === 'OK')?.click();
    }, {description: 'Filtering dataset by substructure', delay: 2000})
    .step('Align by scaffold', async () => {
      filters.close();
      await DG.delay(1000);
      grok.shell.o = tv.dataFrame.col('smiles');
      await DG.delay(2000);
      grok.shell.windows.showHelp = false;
      if (!grok.shell.windows.context.visible)
        grok.shell.windows.context.visible = true;
      await awaitCheck(() => document.getElementsByClassName('grok-entity-prop-panel').length > 0,
        'Property panel did not appear', 5000);
      propPanel = document.getElementsByClassName('grok-entity-prop-panel')[0];
      closeAllAccordionPanes(propPanel!);
      const chemistryPaneContent = getAccordionPane('Chemistry', propPanel!);
      const renderingPaneContent = getAccordionPane('Rendering', chemistryPaneContent!) as HTMLElement;
      await DG.delay(1000);
      const scaffoldSketcher = await openSketcher(renderingPaneContent, 'sketch-link');
      const scaffoldSketcherInput = scaffoldSketcher!
        .getElementsByClassName('grok-sketcher-input')[0]?.children[0] as HTMLInputElement;

      let dT = null;
      try {
        dT = new DataTransfer();
      } catch (e) {}
      const evt = new ClipboardEvent('paste', {clipboardData: dT});
        evt.clipboardData!.setData('text/plain', demoScaffold);
        scaffoldSketcherInput.value = demoScaffold;
        await DG.delay(100);
        scaffoldSketcherInput.dispatchEvent(evt);
        Array.from(scaffoldSketcher!.getElementsByTagName('span')).find((el) => el.textContent === 'OK')?.click();
    }, {description: 'Aligning structures by scaffold', delay: 1000})
    .step('Add sparkline columns', async () => {
      tv.grid.columns.add({gridColumnName: `radar`, cellType: 'radar'});
      tv.grid.columns.add({gridColumnName: `barchart`, cellType: 'barchart'});
      tv.grid.columns.setOrder(firstCols.concat(['radar', 'barchart']).concat(lastCols));
      tv.grid.scrollToCell('MolWt', 0);
    })
    .step('Add color coding', async () => {
      table.col('MolWt')!.meta.colors.setLinear();
      table.col('NOCount')!.meta.colors.setConditional({'0 - 6.25': '#73aff5', '6.25 - 12.50': '#ffa500', '12.50 - 18.75': '#ff5140', '18.75 - 25': '#50af28'});
      table.col('RingCount')!.meta.colors.setConditional();
      grok.shell.windows.showHelp = true;
      //@ts-ignore
      grok.shell.windows.help.showHelp('/help/datagrok/solutions/domains/chem/chem');
      DG.chem.currentSketcherType = sketcherType;
    })
    .start();
}


export async function _demoSimilaritySearch(): Promise<void> {
  const demoScript = new DemoScript('Demo', 'Searching for molecules most similar to target molecule');
  let table: DG.DataFrame;
  let tv: DG.TableView;
  await demoScript
    .step('Load data', async () => {
      tv = await openMoleculeDataset('smiles.csv');
      table = tv.dataFrame;
      grok.shell.windows.showContextPanel = false;
      grok.shell.windows.showHelp = false;
    }, {description: 'Load dataset with molecule columns', delay: 2000})
    .step('Show molecules, most similar to the current', async () => {
      await DG.delay(1000);
      const similarityViewer = tv.addViewer('Chem Similarity Search');
      grok.shell.o = similarityViewer;
    }, {description: 'Open similarity search viewer. Selected molecule becomes target.', delay: 2000})
    .step('Change target molecule', async () => {
      table.currentRowIdx = 2;
      await DG.delay(3000);
      table.currentRowIdx = 10;
      await DG.delay(3000);
      table.currentRowIdx = 3;
    }, {description: 'Fast similarity search re-calculating when changing current molecule', delay: 3000})
    .step('Final', async () => console.log('Finished'))
    .start();
}


export async function _demoSimilarityDiversitySearch(): Promise<void> {
  const tv = await openMoleculeDataset('demo_files/smiles.csv');
  _package.files.readAsText('demo_files/similarity_diversity.layout').then(async (layoutString: string) => {
    const layout = DG.ViewLayout.fromJson(layoutString);
    await DG.delay(100);
    tv.loadLayout(layout);
    grok.shell.windows.showHelp = true;
    setTimeout(() => grok.shell.windows.help.showHelp('/help/datagrok/solutions/domains/chem/chem#similarity-and-diversity-search'), 1000);
  });
}


const SAR_MATRIX_VIEWER = 'SAR Matrix Viewer';

interface SarHint {
  anchor: () => Element | null;
  position: `${ui.hints.POSITION}`;
  title: string;
  text: string;
}

/** A tab header of the SAR viewer, used as the fallback anchor for hints about a tab the tour has
 *  not opened yet — that tab's own markup does not exist until it is shown once. */
function sarTabHeader(title: string): Element | null {
  return Array.from(document.querySelectorAll('.d4-tab-header'))
    .find((e) => e.textContent?.trim() === title) ?? null;
}

/** The matrix is a canvas grid, so its cells are not elements to anchor to. This lays a marker over
 *  the first data cell — past the core column and the substituent header — so the hint about cells
 *  opens against one instead of against the pane. Fixed-positioned and outside the grid so it cannot
 *  disturb the layout it is measuring, and replaced on each call rather than accumulating. */
function matrixCellAnchor(): Element | null {
  const host = document.querySelector('.chem-sar-grid-host');
  if (!(host instanceof HTMLElement))
    return null;
  document.querySelector('.chem-sar-hint-cell')?.remove();
  const box = host.getBoundingClientRect();
  const marker = ui.div([], 'chem-sar-hint-cell');
  marker.style.position = 'fixed';
  marker.style.pointerEvents = 'none';
  marker.style.left = `${box.left + CORE_W}px`;
  marker.style.top = `${box.top + COL_HEADER_H}px`;
  marker.style.width = `${CELL_W}px`;
  marker.style.height = `${CELL_H}px`;
  document.body.appendChild(marker);
  return marker;
}

const SAR_HINTS: SarHint[] = [
  {
    // The first card, not the list: the list runs the full height of the pane, so a popup beside it
    // is centred halfway down and lands over the matrix rather than beside what it describes.
    anchor: () => document.querySelector('.chem-sar-nav-list:not(.chem-sar-xfer-list) .chem-sar-card'),
    position: ui.hints.POSITION.RIGHT,
    title: 'Analog series',
    text: 'Every card on the left is an analog series: compounds that share a core and vary at one ' +
      'position, which is what makes their potencies comparable. The card shows the core, how many ' +
      'cores and compounds it covers, and its score on the current ranking.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-nav-list:not(.chem-sar-xfer-list) .chem-sar-scaffold-card') ??
      document.querySelector('.chem-sar-nav-list:not(.chem-sar-xfer-list) .chem-sar-card'),
    position: ui.hints.POSITION.RIGHT,
    title: 'One scaffold, several series',
    text: 'A Scaffold row gathers the series built on the same scaffold, drawn above them with every ' +
      'position it varies marked. They are listed rather than merged: each varies a different position, ' +
      'and substituents at different positions do not belong in one column.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-nav-controls .ui-input-root'),
    position: ui.hints.POSITION.RIGHT,
    title: 'Rank by',
    text: 'Reorders the list: Potent compounds puts the strongest series first, SAR discontinuity ' +
      'the ones where a small change swings potency most, and Preferred substituent the ones a ' +
      'single R-group dominates.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-nav-controls .chem-sar-struct-icon'),
    position: ui.hints.POSITION.RIGHT,
    title: 'Filter the series',
    text: 'The funnel narrows the series. Filter by Best or Mean potency to keep only the ' +
      'series worth reading, by Spread for the ones with the sharpest SAR, or by Compounds, Cores ' +
      'and Level to drop the thin or over-folded ones. Core searches by substructure.',
  },
  {
    // A card that actually folds: the chevron and the level badge are what the text is about, and a
    // leaf card carries neither. The scaffold row has a chevron but no level, so it is not one either.
    anchor: () => Array.from(document.querySelectorAll(
      '.chem-sar-nav-list:not(.chem-sar-xfer-list) .chem-sar-card:not(.chem-sar-scaffold-card)'))
      .find((card) => card.querySelector('.chem-sar-card-twisty.fa-chevron-down, ' +
        '.chem-sar-card-twisty.fa-chevron-right')) ??
      document.querySelector('.chem-sar-nav-list:not(.chem-sar-xfer-list) .chem-sar-card'),
    position: ui.hints.POSITION.RIGHT,
    title: 'Select and unfold',
    text: 'Click a card to open that series in the matrix on the right. Levels nest: an L1 card is a ' +
      'single core with its substituents, while L2 and L3 group the cores that agree one and two cuts ' +
      'deeper - their badge reads "L2·6" for the L2 matrix containing 6 L1 matrices. Use the chevron to ' +
      'unfold one into the matrices it groups, shown indented beneath it. An orange badge warns that ' +
      'unfolding will not reach every compound: the rest sit in groups too thin to form a matrix of ' +
      'their own, so they are counted on this card but appear in none of the series below.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-chips'),
    position: ui.hints.POSITION.BOTTOM,
    title: 'What the open matrix holds',
    text: 'The chips above the matrix summarise it: how many compounds it holds, its size as cores by ' +
      'substituents, the activity range on the current scale, how many analogs are predicted rather ' +
      'than measured, and how many compounds the dataset holds with no activity value at all. ' +
      'The size counts what is on screen, so a filter shrinks it.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-control-bar .chem-sar-filter-icon'),
    position: ui.hints.POSITION.LEFT,
    title: 'Filter the matrix',
    text: 'The funnel narrows the cells themselves - by Potency, Reference points, MW, Core and substituents R.',
  },
  {
    anchor: () => matrixCellAnchor(),
    position: ui.hints.POSITION.RIGHT,
    title: 'Then click a cell',
    text: 'Cores run down the rows, substituents across the columns. The header over the first column ' +
      'draws the aligned core they all share, and each substituent header carries its R group above ' +
      'the position it occupies, so every column reads against the same core. A solid cell is a compound ' +
      'you have, a dashed one an analog nobody has made. Click any of them to inspect it.',
  },
  {
    anchor: () => document.querySelector('.grok-prop-panel'),
    position: ui.hints.POSITION.LEFT,
    title: 'A compound that exists',
    text: 'Select a solid cell and the Context Panel shows the compound itself: its structure, the ' +
      'potency measured on it, and the compounds sharing its core and its substituent as the SAR ' +
      'around it. Some solid cells carry a "~" value instead of a plain one - those are compounds ' +
      'the dataset already holds but dont have any activity. They keep the solid frame and their ' +
      'registration id, because the compound is real and only the number is estimated.',
  },
  {
    anchor: () => document.querySelector('.grok-prop-panel'),
    position: ui.hints.POSITION.LEFT,
    title: 'A virtual analog',
    text: 'Select a dashed cell and the Context Panel explains it: the core and substituent were ' +
      'never combined, so its potency is predicted from how each performs elsewhere in the matrix, ' +
      'with the compounds behind every term listed as reference points. Add to Make list collects ' +
      'the ones worth making into the Make list tab.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-control-bar .chem-sar-cart-icon'),
    position: ui.hints.POSITION.LEFT,
    title: 'Add to the Make list',
    text: 'The cart adds whichever cell is selected in the matrix to the Make list, so an analog can ' +
      'be collected without leaving the grid. It takes measured compounds as well as predicted ones.',
  },
  {
    anchor: () => sarTabHeader('SAR Transfer'),
    position: ui.hints.POSITION.BOTTOM,
    title: 'Carry the SAR across scaffolds',
    text: 'The SAR Transfer tab pairs cores whose potency trends run in parallel over the substituents ' +
      'they share. Where the correlation holds, an optimization found on one scaffold is expected to ' +
      'carry to the other, and the analogs it argues for are marked in the matrix.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-xfer-nav .chem-sar-scaffold-card') ??
      document.querySelector('.chem-sar-xfer-list') ?? sarTabHeader('SAR Transfer'),
    position: ui.hints.POSITION.RIGHT,
    title: 'Sources, nested by series',
    text: 'Each card on the left is one source core, and the cards are gathered under the series that ' +
      'core belongs to. Use the chevron to fold a series shut while you read another. A source that transfers ' +
      'to several targets is one card rather than several: the dropdown in the header above the grid ' +
      'switches between its targets, keeping the source fixed.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-xfer-nav .chem-sar-nav-controls .ui-input-root') ??
      sarTabHeader('SAR Transfer'),
    position: ui.hints.POSITION.RIGHT,
    title: 'Rank the transfers',
    text: 'Potent compounds orders by the strongest analog a pairing points at - the question the ' +
      'matrix list asks, asked of a transfer. Fold match orders by how far the size of each step ' +
      'carries rather than merely its direction: two cores can rank every substituent alike and still ' +
      'disagree on what a swap is worth. Analogs gained orders by how many compounds the pairing ' +
      'actually argues for making, which is none when both cores have already explored the same ' +
      'R-groups. The funnel beside the dropdown filters on all three.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-xfer-panel .chem-sar-grid-host') ??
      sarTabHeader('SAR Transfer'),
    position: ui.hints.POSITION.LEFT,
    title: 'Analogs the transfer argues for',
    text: 'The two series are laid out over the substituents they share, and the columns captioned ' +
      '"predicted" are the ones only one of them has tried. Those cells are dashed, exactly as in the ' +
      'matrix: click one and the Context Panel gives the same Free-Wilson breakdown and the same Add ' +
      'to Make list, so an analog suggested by a transfer joins the same Make list.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-xfer-panel .chem-sar-control-bar .ui-input-root') ??
      document.querySelector('.chem-sar-xfer-panel .chem-sar-main-bar') ?? sarTabHeader('SAR Transfer'),
    position: ui.hints.POSITION.BOTTOM,
    title: 'Pick the target core',
    text: 'One core often transfers to several. The dropdown in the header switches between them ' +
      'while the source stays put, so you can read the same series against each core its SAR carries ' +
      'to in turn. The title left of it always names the pair on screen, and the chips beside it - ' +
      'correlation, fold match, shared R-groups and how many compounds the pairing argues for - ' +
      'change with the target you pick. A source with only one target gets no dropdown.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-xfer-panel .chem-sar-cart-icon') ??
      document.querySelector('.chem-sar-xfer-panel .chem-sar-main-bar') ?? sarTabHeader('SAR Transfer'),
    position: ui.hints.POSITION.LEFT,
    title: 'Collect from a transfer',
    text: 'The cart adds whichever cell is selected to the Make list, the same as on the matrix tab. ' +
      'Click a predicted cell on either side to take the analog that pairing argues for.',
  },
  {
    anchor: () => sarTabHeader('Make list'),
    position: ui.hints.POSITION.BOTTOM,
    title: 'The Make list',
    text: 'Everything collected lands here: each analog with its structure, its potency and the ' +
      'activity that potency is read on, how much evidence stood behind it, and the series, core and ' +
      'substituent it came from. Status and Method are separate columns because they answer different ' +
      'questions - Status says whether the compound exists (Synthesized, Untested or Virtual), Method ' +
      'says where its number came from (measured or predicted), and an untested compound is made and ' +
      'predicted at once. Add to workspace hands out a copy to save, export or join; Remove drops the ' +
      'selected compound and Clear empties the list.',
  },
];

function showSarHint(i: number, previous?: HTMLElement): void {
  previous?.remove();
  if (i >= SAR_HINTS.length) {
    document.querySelector('.chem-sar-hint-cell')?.remove();
    return;
  }
  const hint = SAR_HINTS[i];
  const anchor = hint.anchor();
  if (!(anchor instanceof HTMLElement)) {
    showSarHint(i + 1);
    return;
  }
  // Cards sit wherever their rank puts them, and one low in the list leaves too little height beside
  // it for the popup, which then runs past the bottom of the pane. Lifting it to the top of the
  // scroller first gives the popup the full pane to open into.
  if (anchor.closest('.chem-sar-nav-list'))
    anchor.scrollIntoView({block: 'start'});
  const content = ui.div();
  content.append(ui.h3(hint.title));
  content.append(ui.divText(hint.text));
  const buttonHost = ui.divH([], {style: {justifyContent: 'flex-end'}});
  content.append(buttonHost);
  const popup = ui.hints.addHint(anchor, content, hint.position);
  buttonHost.append(i === SAR_HINTS.length - 1 ?
    ui.button('Close', () => {
      popup.remove();
      document.querySelector('.chem-sar-hint-cell')?.remove();
    }) :
    ui.button('Next', () => showSarHint(i + 1, popup)));
}

function showSarHints(): void {
  const timer = setInterval(() => {
    if (document.querySelector('.chem-sar-grid-host')) {
      clearInterval(timer);
      showSarHint(0);
    }
  }, 250);
  setTimeout(() => clearInterval(timer), 120000);
}

export async function _demoSarMatrix(): Promise<void> {
  grok.shell.windows.showContextPanel = true;
  grok.shell.windows.showHelp = true;
  const p = await grok.functions.eval('Chem:SarMatrixDemo');
  const project = await grok.dapi.projects.find(p.id);
  await project.open();
  await DG.delay(300);
  const tv = grok.shell.tv;
  if (tv) {
    for (const viewer of tv.viewers) {
      if (viewer.type === SAR_MATRIX_VIEWER)
        dockSarMatrixTabs(tv, viewer);
    }
  }
  showSarHints();
  setTimeout(() => grok.shell.windows.help.showHelp('/help/datagrok/solutions/domains/chem/chem#sar-matrix'), 1000);
}

export async function _demoMMPA(): Promise<void> {
  const tv = await openMoleculeDataset('demo_files/mmp_demo.csv');

  _package.files.readAsText('demo_files/mmp_demo.layout').then(async (layoutString: string) => {
    const layout = DG.ViewLayout.fromJson(layoutString);
    tv.loadLayout(layout);
    tv.dataFrame.currentRowIdx = 0;
    try {
      await awaitCheck(() => {
        for (const v of tv.viewers) {
          if (v.type === 'Matched Molecular Pairs Analysis') {
            v.helpUrl = 'https://raw.githubusercontent.com/datagrok-ai/public/refs/heads/master/help/datagrok/solutions/domains/chem/chem.md#matched-molecular-pairs';
            return true;
          }
        }
        return false;
      }, '', 20000);
    } catch (e) {};
    setTimeout(()=> {
      grok.shell.windows.showHelp = true;
      grok.shell.windows.help.showHelp('/help/datagrok/solutions/domains/chem/chem#matched-molecular-pairs');
    }, 1000);

    // grok.shell.windows.showHelp = true;
    // grok.shell.windows.help.showHelp('/help/datagrok/solutions/domains/chem/#matched-molecular-pairs');
  });
}


export async function _demoRgroupAnalysis(): Promise<void> {
  const demoScript = new DemoScript('R-Group Analysis', 'Performing R Group Analysis',
    undefined, {autoStartFirstStep: true});
  let table: DG.DataFrame;
  let tv: DG.TableView;
  let sketcherInput: HTMLInputElement;
  let sketcher: Element;

  const findTrellisPlot = () => {
    for (const viewer of tv.viewers) {
      if (viewer.type === DG.VIEWER.TRELLIS_PLOT)
        return viewer;
    }
    return null;
  };

  await demoScript
    .step('Load data', async () => {
      tv = await openMoleculeDataset('demo_files/sar_small.csv');
      table = tv.dataFrame;
      grok.shell.windows.showContextPanel = false;
      grok.shell.windows.showHelp = false;
    }, {description: 'Load dataset with molecule columns', delay: 2000})
    .step('Specify scaffold', async () => {
      await DG.delay(1000);
      rGroupAnalysis(table.col('smiles')!);
      await DG.delay(3000);
      sketcher = document.getElementsByClassName('d4-dialog')[0];
      sketcherInput = sketcher!.getElementsByClassName('grok-sketcher-input')[0]?.children[0] as HTMLInputElement;
      sketcherInput.value = 'O=C1CN=C(c2ccccc2N1)C3CCCCC3';
      const progressBar = DG.TaskBarProgressIndicator.create(`Sketcher initialization in progress...`);
      await DG.delay(3000);
      progressBar.close();
      sketcherInput.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter'}));
    }, {description: 'Open R Group Analysis viewer and enter scaffold structure', delay: 2000})
    .step('Analyse R Groups', async () => {
      const dlgOKButton = Array.from(sketcher!.getElementsByTagName('span')).find((el) => el.textContent === 'OK');
      if (dlgOKButton)
        dlgOKButton.click();
      await awaitCheck(() => {
        return !!findTrellisPlot();
      },
      'r group analysis has not been loaded', 30000);
    }, {description: 'Trellis plot is created from R Group Analysis results', delay: 2000})
    .step('Explore results in various viewers', async () => {
      await DG.delay(1000);
      tv.scatterPlot({x: 'R1', y: 'R2', jitterSize: 4, size: 'LD(50)', color: 'Mol Wt.', autoAxisSize: false});
      tv.barChart({split: 'R1'});
    }, {description: 'Any other type of viewer can be easily created on R Group analysis results', delay: 2000})
    .start();
}


export async function _demoActivityCliffs(): Promise<void> {
  const demoScript = new DemoScript('Activity Cliffs',
    'Searching similar structures with significant activity difference', undefined, {autoStartFirstStep: true});
  let table: DG.DataFrame;
  let tv: DG.TableView;
  let scatterPlot: DG.Viewer;
  await demoScript
    .step('Load data', async () => {
      tv = await openMoleculeDataset('demo_files/sar_small.csv');
      table = tv.dataFrame;
    }, {description: 'Load dataset with molecule and activity columns.', delay: 2000})
    .step('Find activity cliffs', async () => {
      const molecules = table.col('smiles')!;

      const preprocessing = DG.Func.find({name: 'getFingerprints', package: 'Chem'})[0];
      await PackageFunctions.activityCliffs(table, molecules, table.col('In-vivo Activity')!,
        78, DimReductionMethods.T_SNE, BitArrayMetricsNames.Tanimoto,
        preprocessing, {}, true);
      // tv = (grok.shell.view('Browse')! as DG.BrowseView)!.preview! as DG.TableView;
      awaitCheck(() => {
        for (const v of tv.viewers) {
          if (v.type === DG.VIEWER.SCATTER_PLOT) {
            scatterPlot = v;
            return true;
          }
        }
        return false;
      }, '', 10000);
      await DG.delay(1000);
    }, {description: `Results are shown on a scatter plot. Each point on a scatter plot corresponds to a molecule from a dataset.
    Pairs of molecules with similarity higher than specified cutoff, are connected by lines. Marker color corresponds to molecule activity.
    Line opacity corresponds to molecule pair SALI value (Structure−Activity Landscape Index - activity difference divided by 1 minus similarity).
    Marker size corresponds to highest SALI value detected for the molecule.`, delay: 2000})
    .step('Explore activity cliffs', async () => {
      await DG.delay(1000);
      const link = Array.from(scatterPlot!.root.getElementsByClassName('scatter_plot_link'));
      if (link.length)
        (link[0] as HTMLElement).click();
      await DG.delay(1000);
    }, {description: `Detected cliffs are available in a separate table. 
    Cliffs are pairs of molecules with similarity higher than cutoff. Cliffs are sorted by SALI value.`, delay: 2000})
    .step('Select cliffs', async () => {
      await DG.delay(1000);
      let cliffsGrid: DG.Viewer | null = null;
      for (const i of tv.viewers) {
        if (i.dataFrame.name === `${CLIFFS_DF_NAME}${activityCliffsIdx}`)
          cliffsGrid = i;
      }
      cliffsGrid!.dataFrame.currentRowIdx = 35;
      await DG.delay(3000);
      cliffsGrid!.dataFrame.currentRowIdx = 6;
      await DG.delay(3000);
      cliffsGrid!.dataFrame.currentRowIdx = 5;
    }, {description: `To zoom scatter plot to exact cliff, click on a row in the cliffs table. 
    Additional information about molecule pair is on the context panel. Non common fragments are highlighted in molecules.`, delay: 3000})
    .start();
}

export async function _demoActivityCliffsLayout(): Promise<void> {
  grok.shell.windows.showContextPanel = true;
  grok.shell.windows.showHelp = true;
  const p = await grok.functions.eval('Chem:DemoActivityCliffs');
  const project = await grok.dapi.projects.find(p.id);
  await project.open();
  let scatterPlot: DG.Viewer | null = null;
  for (const i of grok.shell.tv.viewers) {
    if (i.type == DG.VIEWER.SCATTER_PLOT)
      scatterPlot = i;
  }
  let cliffsLink;
  try {
    await awaitCheck(() => {
      const link = scatterPlot?.root.getElementsByClassName('scatter_plot_link');
      if (link?.length) {
        cliffsLink = link[0];
        return true;
      }
      return false;
    }, '', 10000);
    (cliffsLink as any as HTMLElement).click();
  } catch (e) {}
  setTimeout(() => grok.shell.windows.help.showHelp('/help/datagrok/solutions/domains/chem/chem#activity-cliffs'), 1000);
}

export async function _demoRGroups(): Promise<void> {
  grok.shell.windows.showContextPanel = true;
  grok.shell.windows.showHelp = true;
  const p = await grok.functions.eval('Chem:RGroupsDemo');
  const project = await grok.dapi.projects.find(p.id);
  await project.open();
  setTimeout(() => grok.shell.windows.help.showHelp('/help/datagrok/solutions/domains/chem/chem#r-groups-analysis'), 1000);
}

export async function _demoChemicalSpace(): Promise<void> {
  grok.shell.windows.showContextPanel = true;
  grok.shell.windows.showHelp = true;
  const p = await grok.functions.eval('Chem:ChemicalSpaceDemo');
  const project = await grok.dapi.projects.find(p.id);
  await project.open();
  grok.functions.call('Dendrogram:HierarchicalClustering', {
    df: grok.shell.project.children.find((it) => it instanceof DG.TableInfo)?.dataFrame,
    colNameList: ['molecule'],
    distance: 'euclidian',
    linkage: 'average',
  });
  //for dendrogram to render correctly
  const sub = grok.shell.tv?.grid?.onAfterDrawOverlay?.subscribe(() => {
    sub.unsubscribe();
    setTimeout(() => grok.shell.tv?.grid?.invalidate(), 100);
  });
  setTimeout(() => grok.shell.windows.help.showHelp('/help/datagrok/solutions/domains/chem/chem#chemical-space'), 1000);
}

export async function _demoScaffoldTree(): Promise<void> {
  const tv = await openMoleculeDataset('mol1K.csv');
  _package.files.readAsText('demo_files/mol1K.layout').then(async (layoutString: string) => {
    const layout = DG.ViewLayout.fromJson(layoutString);
    tv.loadLayout(layout);

    const {dataFrame} = tv;
    const scaffoldTree = (await dataFrame.plot.fromType(DG.VIEWER.SCAFFOLD_TREE)) as unknown as ScaffoldTreeViewer;
    scaffoldTree.root.classList.add('d4-viewer', 'd4-scaffold-tree');
    const treeStr = await _package.files.readAsText('demo_files/scaffold-tree.json');
    const table: DG.DataFrame = tv.dataFrame;
    await grok.data.detectSemanticTypes(table);

    scaffoldTree.molCol = table.columns.bySemType(DG.SEMTYPE.MOLECULE);
    scaffoldTree.dataFrame = table;
    scaffoldTree.size = 'small';

    tv.dockManager.dock(scaffoldTree, DG.DOCK_TYPE.LEFT, null, 'Scaffold Tree', 0.4);
    await scaffoldTree.loadTreeStr(treeStr);

    grok.shell.windows.showHelp = true;
    setTimeout(() => grok.shell.windows.help.showHelp('/help/datagrok/solutions/domains/chem/chem#scaffold-tree-analysis'), 1000);
  });
}
