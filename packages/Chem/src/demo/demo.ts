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
  await _package.files.readAsText('demo_files/similarity_diversity.layout').then(async (layoutString: string) => {
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
  position?: `${ui.hints.POSITION}`;
  title: string;
  text: string;
  /** Put the viewer where the hint's anchor exists — a tab, a section — before it is looked for. The
   *  tour walks four tabs and three sections of one of them, and an anchor on a pane that is not open
   *  is not in the document at all. */
  before?: () => void;
}

/** The one element matching `selector` whose text is exactly `title`. */
function byText(selector: string, title: string): HTMLElement | null {
  const found = Array.from(document.querySelectorAll(selector))
    .find((e) => e.textContent?.trim() === title);
  return found instanceof HTMLElement ? found : null;
}

/** Click a tab header of the SAR viewer, or of the Summary tab's own sections, by name. */
function sarTab(title: string): void {
  sarTabHeader(title)?.click();
}

/** A tab header of the SAR viewer, used as the fallback anchor for hints about a tab the tour has
 *  not opened yet — that tab's own markup does not exist until it is shown once. */
function sarTabHeader(title: string): HTMLElement | null {
  return byText('.d4-tab-header', title);
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

/** Where the two tours point. Both walk the same analysis, so a selector belongs in one place: the
 *  fallback is what the step lands on while its pane is still building. */
const at = {
  scale: () => document.querySelector('.chem-sar-sum-scale') ?? sarTabHeader('Summary'),
  coverage: () => document.querySelector('.chem-sar-sum-coverage') ??
    document.querySelector('.chem-sar-sum-overview'),
  components: () => document.querySelector('.chem-sar-sum-comp') ?? document.querySelector('.chem-sar-sum-answers'),
  swapRow: () => document.querySelector('.chem-sar-sum-swap') ?? document.querySelector('.chem-sar-sum-answers'),
  making: () => document.querySelector('.chem-sar-sum-making-block') ?? document.querySelector('.chem-sar-sum-making'),
  navCard: () => document.querySelector('.chem-sar-nav-list:not(.chem-sar-xfer-list) .chem-sar-card') ??
    sarTabHeader('SAR Matrix'),
  chips: () => document.querySelector('.chem-sar-chips'),
  contextPanel: () => document.querySelector('.grok-prop-panel'),
  transferTab: () => sarTabHeader('SAR Transfer'),
  transferSources: () => document.querySelector('.chem-sar-xfer-nav .chem-sar-scaffold-card') ??
    document.querySelector('.chem-sar-xfer-list') ?? sarTabHeader('SAR Transfer'),
  transferGrid: () => document.querySelector('.chem-sar-xfer-panel .chem-sar-grid-host') ??
    document.querySelector('.chem-sar-xfer-panel') ?? sarTabHeader('SAR Transfer'),
  makeList: () => sarTabHeader('Make list'),
};

const SAR_HINTS: SarHint[] = [
  {
    before: () => sarTab('Summary'),
    anchor: at.scale,
    title: 'Start on the Summary',
    text: 'What the analysis concluded. This line is the scale everything rests on: the endpoint, its ' +
      'transform, which direction is better and the observed range.',
  },
  {
    anchor: at.coverage,
    title: 'What reached a matrix',
    text: 'How much of the table reached a series: measured, never assayed, or assayed with no analog ' +
      'to pair with.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-sum-list .chem-sar-card') ??
      document.querySelector('.chem-sar-sum-overview'),
    position: ui.hints.POSITION.TOP,
    title: 'Start here',
    text: 'The series to open first. The squares say why: most compounds, widest range, best ' +
      'compound, most predictions, best-validated fit. Click a row for detail; its arrow opens the matrix.',
  },
  {
    before: () => at.components()?.scrollIntoView({block: 'nearest'}),
    anchor: at.components,
    title: 'What to change',
    text: 'One additive (Free-Wilson) fit over the whole table. Each component shows its spread, its ' +
      'best value and how far it sits above the average. The widest spread is the one to change next.',
  },
  {
    anchor: at.swapRow,
    title: 'Best measured swap, per component',
    text: 'Measured pairs that differ in one component only, with no model involved. Click one for its ' +
      'full pool.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-sum-answers') ?? document.querySelector('.chem-sar-sum-comp'),
    title: 'The evidence, in place',
    text: 'Click a tile to open the evidence behind its answer right underneath; click it again to close it.',
  },
  {
    before: () => document.querySelector('.chem-sar-sum-making')?.scrollIntoView(),
    anchor: at.making,
    position: ui.hints.POSITION.TOP,
    title: 'Worth making',
    text: 'Click to open: analogs not in the table, ranked on their gain over the best compound their ' +
      'series made. Compounds you hold but never assayed sit just above.',
  },
  {
    before: () => sarTab('SAR Matrix'),
    anchor: at.navCard,
    position: ui.hints.POSITION.RIGHT,
    title: 'The series themselves',
    text: 'The evidence behind the Summary: the series, their matrices, the transfers and the Make ' +
      'list.',
  },
  {
    // The first card, not the list: the list runs the full height of the pane, so a popup beside it
    // is centred halfway down and lands over the matrix rather than beside what it describes.
    anchor: () => document.querySelector('.chem-sar-nav-list:not(.chem-sar-xfer-list) .chem-sar-card'),
    position: ui.hints.POSITION.RIGHT,
    title: 'Analog series',
    text: 'Each card is an analog series: compounds that share a core and vary at one position. It ' +
      'shows the core, its counts and its score.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-nav-list:not(.chem-sar-xfer-list) .chem-sar-scaffold-card') ??
      document.querySelector('.chem-sar-nav-list:not(.chem-sar-xfer-list) .chem-sar-card'),
    position: ui.hints.POSITION.RIGHT,
    title: 'One scaffold, several series',
    text: 'A Scaffold row groups the series built on one scaffold. They stay separate because each ' +
      'varies a different position.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-nav-controls .ui-input-root'),
    position: ui.hints.POSITION.RIGHT,
    title: 'Rank by',
    text: 'Reorders the list by potency, by SAR discontinuity, or by a dominant substituent.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-nav-controls .chem-sar-struct-icon'),
    position: ui.hints.POSITION.RIGHT,
    title: 'Filter the series',
    text: 'Narrows the series by potency, spread, size, level or core substructure.',
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
    text: 'Click a card to open its matrix. L2 and L3 group cores that agree one or two cuts deeper ' +
      '("L2·6" holds six L1 matrices); the chevron unfolds them. An orange badge means some ' +
      'compounds sit in groups too thin to show.',
  },
  {
    anchor: at.chips,
    title: 'What the open matrix holds',
    text: 'Compounds, size as cores × substituents, activity range, predicted analogs and compounds ' +
      'with no value. The size follows any filter.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-control-bar .chem-sar-filter-icon'),
    position: ui.hints.POSITION.LEFT,
    title: 'Filter the matrix',
    text: 'Narrows the cells by potency, reference points, MW, core or substituent.',
  },
  {
    anchor: matrixCellAnchor,
    position: ui.hints.POSITION.RIGHT,
    title: 'Then click a cell',
    text: 'Cores run down the rows, substituents across. A solid cell is a compound you have, a ' +
      'dashed one an analog nobody has made. Click one to inspect it.',
  },
  {
    anchor: at.contextPanel,
    position: ui.hints.POSITION.LEFT,
    title: 'A compound that exists',
    text: 'Its structure, measured potency and SAR neighbours. A "~" value marks a compound you hold ' +
      'without a measurement: the compound is real, the number estimated.',
  },
  {
    anchor: at.contextPanel,
    position: ui.hints.POSITION.LEFT,
    title: 'A virtual analog',
    text: 'Never made, so its potency is predicted from how its core and substituent do elsewhere, ' +
      'with the compounds behind each term. Add to Make list collects it.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-control-bar .chem-sar-cart-icon'),
    position: ui.hints.POSITION.LEFT,
    title: 'Add to the Make list',
    text: 'Adds the selected cell, measured or predicted, to the Make list.',
  },
  {
    anchor: at.transferTab,
    title: 'Carry the SAR across scaffolds',
    text: 'Pairs cores whose potency trends run in parallel, so an optimization found on one may ' +
      'carry to the other.',
  },
  {
    anchor: at.transferSources,
    position: ui.hints.POSITION.RIGHT,
    title: 'Sources, nested by series',
    text: 'One card per source core, grouped by series. A source with several targets switches ' +
      'between them in the dropdown above the grid.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-xfer-nav .chem-sar-nav-controls .ui-input-root') ??
      sarTabHeader('SAR Transfer'),
    position: ui.hints.POSITION.RIGHT,
    title: 'Rank the transfers',
    text: 'Order by Potent compounds, Fold match (whether step sizes carry, not just direction) or ' +
      'Analogs gained. The funnel filters on all three.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-xfer-panel .chem-sar-grid-host') ??
      sarTabHeader('SAR Transfer'),
    position: ui.hints.POSITION.LEFT,
    title: 'Analogs the transfer argues for',
    text: 'Both series over their shared substituents. Dashed "predicted" cells are what the transfer ' +
      'suggests; click one for its breakdown and Add to Make list.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-xfer-panel .chem-sar-control-bar .ui-input-root') ??
      document.querySelector('.chem-sar-xfer-panel .chem-sar-main-bar') ?? sarTabHeader('SAR Transfer'),
    title: 'Pick the target core',
    text: 'Switches the target core while the source stays put; the chips follow the target.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-xfer-panel .chem-sar-cart-icon') ??
      document.querySelector('.chem-sar-xfer-panel .chem-sar-main-bar') ?? sarTabHeader('SAR Transfer'),
    position: ui.hints.POSITION.LEFT,
    title: 'Collect from a transfer',
    text: 'Adds the selected cell to the Make list, as on the matrix tab.',
  },
  {
    anchor: at.makeList,
    title: 'The Make list',
    text: 'Everything you collected, with structure, potency, evidence and origin. Status says ' +
      'whether it exists, Method whether its number is measured or predicted.',
  },
];

/**
 * The PROTAC tour. Shorter than the fragmentation one and about a different thing: here the three
 * parts of every molecule are already columns, so nothing is cut and the question is which of the
 * three to change.
 */
const PROTAC_HINTS: SarHint[] = [
  {
    before: () => sarTab('Summary'),
    anchor: at.scale,
    title: '2 792 degraders, split three ways',
    text: 'Warhead, linker and E3 ligand are columns of their own. The endpoint is predicted ' +
      'solubility in log units, higher is better.',
  },
  {
    anchor: at.coverage,
    title: 'What reached a matrix',
    text: 'How many degraders landed in a series, and how many had no analog to pair with.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-sum-list'),
    position: ui.hints.POSITION.TOP,
    title: 'Which series to open first',
    text: 'The series to open first; the squares say why each is listed. Click a row for detail; its ' +
      'arrow opens the matrix.',
  },
  {
    before: () => at.components()?.scrollIntoView({block: 'nearest'}),
    anchor: at.components,
    title: 'Which of the three to change',
    text: 'One additive fit over all three components ranks them by how far each moves solubility, ' +
      'each adjusted for the other two. On the right: the linker is the core, warheads run across the ' +
      'columns and the E3 ligand folds into the row.',
  },
  {
    anchor: at.swapRow,
    title: 'And what was actually measured',
    text: 'Measured pairs alike in every part but one, per component. No model involved.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-sum-comp .chem-sar-sum-half') ?? at.components(),
    title: 'Every linker, warhead and ligand ranked',
    text: 'Click a component\'s value for its full ranking, drawn as structures, or its swap for the measured ' +
      'pairs behind it. Both open right under the row.',
  },
  {
    before: () => document.querySelector('.chem-sar-sum-making')?.scrollIntoView(),
    anchor: at.making,
    position: ui.hints.POSITION.TOP,
    title: 'Combinations nobody has built',
    text: 'Click to open: warhead, linker and E3 ligand combinations never made, ranked on gain over ' +
      'their series\' best. Listed even when none clears the trust gate.',
  },
  {
    before: () => sarTab('SAR Matrix'),
    anchor: at.navCard,
    position: ui.hints.POSITION.RIGHT,
    title: 'One linker per series',
    text: 'Each card is one linker: warheads across the columns, the E3 ligand folded into the row. ' +
      'Best and Mean on the right are its solubility.',
  },
  {
    anchor: at.chips,
    title: 'What the open matrix holds',
    text: 'Degraders measured, cells the components could fill, predictions and the series\' R². ' +
      'Columns, on the right, puts another component across the top.',
  },
  {
    anchor: matrixCellAnchor,
    position: ui.hints.POSITION.RIGHT,
    title: 'Click a cell',
    text: 'One warhead against one E3 ligand on this linker, coloured by solubility.',
  },
  {
    anchor: at.contextPanel,
    position: ui.hints.POSITION.LEFT,
    title: 'A degrader that exists',
    text: 'The whole molecule assembled from its three parts, with its measured value. Nothing here ' +
      'is modelled.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-control-bar .chem-sar-cart-icon') ??
      document.querySelector('.chem-sar-control-bar'),
    position: ui.hints.POSITION.LEFT,
    title: 'Collect what is worth making',
    text: 'A faint cell is a degrader nobody has built. The cart adds the selected cells to the Make ' +
      'list with their predicted values.',
  },
  {
    before: () => sarTab('SAR Transfer'),
    anchor: at.transferTab,
    title: 'Carry the SAR across linkers',
    text: 'Linkers tried with the same warheads whose solubility moves together inform each other.',
  },
  {
    anchor: at.transferSources,
    position: ui.hints.POSITION.RIGHT,
    title: 'Sources, nested by series',
    text: 'Each group is a linker lending its trend. "r 1.00 over 10 shared": ten warheads ranked the ' +
      'same way on both.',
  },
  {
    anchor: at.transferGrid,
    position: ui.hints.POSITION.LEFT,
    title: 'What the transfer argues for',
    text: 'Both matrices side by side. A source value over a target hole is an analog the transfer ' +
      'suggests.',
  },
  {
    before: () => sarTab('Make list'),
    anchor: at.makeList,
    title: 'The Make list',
    text: 'Everything collected, with its three parts, predicted value and origin. It exports like ' +
      'any table.',
  },
];

/** The two SAR demos share this walker and differ only in the step list they hand it. */
function showSarHint(hints: SarHint[], i: number, previous?: HTMLElement, tries = 0): void {
  previous?.remove();
  if (i >= hints.length) {
    document.querySelector('.chem-sar-hint-cell')?.remove();
    return;
  }
  const hint = hints[i];
  if (tries === 0)
    hint.before?.();
  const anchor = hint.anchor();
  if (!(anchor instanceof HTMLElement)) {
    // A pane the step just opened builds on a timer, so its anchor is not in the document the instant
    // the tab is clicked. Waiting a few beats for it is the difference between a step and a skipped
    // step; past that the anchor is genuinely absent and the tour moves on.
    if (tries < 8) {
      setTimeout(() => showSarHint(hints, i, undefined, tries + 1), 250);
      return;
    }
    showSarHint(hints, i + 1);
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
  const popup = ui.hints.addHint(anchor, content, hint.position ?? ui.hints.POSITION.BOTTOM);
  buttonHost.append(i === hints.length - 1 ?
    ui.button('Close', () => {
      popup.remove();
      document.querySelector('.chem-sar-hint-cell')?.remove();
    }) :
    ui.button('Next', () => showSarHint(hints, i + 1, popup)));
}

function showSarHints(hints: SarHint[]): void {
  const timer = setInterval(() => {
    // Either landing surface, not the matrix grid alone: the viewer opens on the Summary, and that
    // grid does not exist until the SAR Matrix tab is shown — which is a step of the tour, not a
    // precondition for starting it.
    if (document.querySelector('.chem-sar-sum-overview') || document.querySelector('.chem-sar-grid-host')) {
      clearInterval(timer);
      showSarHint(hints, 0);
    }
  }, 250);
  setTimeout(() => clearInterval(timer), 120000);
}

/** Restored from a project rather than run: decomposing the table takes long enough to watch a loader. */
async function openSarDemo(projectFunc: string, hints: SarHint[]): Promise<void> {
  grok.shell.windows.showContextPanel = true;
  const p = await grok.functions.eval(projectFunc);
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
  showSarHints(hints);
  setTimeout(() => {
    grok.shell.windows.showHelp = true;
    grok.shell.windows.help.showHelp('/help/datagrok/solutions/domains/chem/chem#sar-matrix');
  }, 1000);
}

export async function _demoSarMatrix(): Promise<void> {
  await openSarDemo('Chem:SarMatrixDemo', SAR_HINTS);
}

/**
 * The same analysis over components that are already columns: 2 792 PROTAC degraders from patents,
 * each split into warhead, linker and E3 ligand. Nothing is cut: the matrices are built by grouping,
 * every component is ranked in one fit, and every component's measured pairs are pooled at once.
 */
export async function _demoSarMatrixProtac(): Promise<void> {
  await openSarDemo('Chem:SarMatrixProtacDemo', PROTAC_HINTS);
}

export async function _demoMMPA(): Promise<void> {
  const tv = await openMoleculeDataset('demo_files/mmp_demo.csv');

  await _package.files.readAsText('demo_files/mmp_demo.layout').then(async (layoutString: string) => {
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
  await _package.files.readAsText('demo_files/mol1K.layout').then(async (layoutString: string) => {
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
