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

const SAR_HINTS: SarHint[] = [
  {
    anchor: () => document.querySelector('.chem-sar-nav-list'),
    position: ui.hints.POSITION.RIGHT,
    title: 'Analog series',
    text: 'Every card on the left is an analog series: compounds that share a core and vary at one ' +
      'position, which is what makes their potencies comparable. The card shows the core, how many ' +
      'cores and compounds it covers, and its score on the current ranking.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-scaffold-card') ??
      document.querySelector('.chem-sar-nav-list'),
    position: ui.hints.POSITION.RIGHT,
    title: 'One scaffold, several series',
    text: 'A Scaffold row gathers the series built on the same scaffold, drawn above them with every ' +
      'position it varies marked. They are listed rather than merged: each varies a different position, ' +
      'and substituents at different positions do not belong in one column. Open them in turn to read ' +
      'the same compounds from each position in turn.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-nav-controls .ui-input-root'),
    position: ui.hints.POSITION.RIGHT,
    title: 'Rank by',
    text: 'Reorders the list: Potent compounds puts the strongest series first, SAR discontinuity ' +
      'the ones where a small change swings potency most, and Preferred substituent the ones a ' +
      'single R-group dominates. The score on each card follows the scheme - best and mean potency ' +
      'for Potent compounds, best R for Preferred substituent.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-nav-controls .chem-sar-struct-icon'),
    position: ui.hints.POSITION.RIGHT,
    title: 'Filter the series',
    text: 'This funnel narrows the list on the left. Filter by Best or Mean potency to keep only the ' +
      'series worth reading, by Spread for the ones with the sharpest SAR, or by Compounds, Cores ' +
      'and Level to drop the thin or over-folded ones. Core searches by substructure.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-nav-list'),
    position: ui.hints.POSITION.RIGHT,
    title: 'Select and unfold',
    text: 'Click a card to open that series in the matrix on the right. Levels nest: an L1 card is a ' +
      'single core with its substituents, while L2 and L3 group the cores that agree one and two cuts ' +
      'deeper - their badge reads "L2·6" for a level-2 matrix folding six in. Use the chevron to ' +
      'unfold one into the matrices it groups, shown indented beneath it. An orange badge warns that ' +
      'unfolding will not reach every compound: the rest sit in groups too thin to form a matrix of ' +
      'their own, so they are counted on this card but appear in none of the series below. The chip ' +
      'above the matrix repeats the split, and both give the exact numbers on hover.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-chips'),
    position: ui.hints.POSITION.BOTTOM,
    title: 'What the open matrix holds',
    text: 'The chips above the matrix summarise it: how many compounds it holds, its size as cores by ' +
      'substituents, the activity range on the current scale, and how many analogs are predicted ' +
      'rather than measured. The size counts what is on screen, so a filter shrinks it.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-control-bar .chem-sar-filter-icon'),
    position: ui.hints.POSITION.LEFT,
    title: 'Filter the matrix',
    text: 'This funnel filters the cells themselves. Reference points is the one to reach for - it ' +
      'hides predictions resting on too little evidence, so what stays is what the data supports. ' +
      'Potency, MW, Core and R narrow the grid the same way, and the Caption dropdown in the same bar ' +
      'annotates each substituent column with a metric - its mean potency or molecular weight.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-grid-host'),
    position: ui.hints.POSITION.RIGHT,
    title: 'Then click a cell',
    text: 'Cores run down the rows, substituents across the columns. The header over the first column ' +
      'draws the aligned core they all share, and each substituent header carries its R group above ' +
      'the position it occupies, so every column reads against the same core. A solid cell is a ' +
      'measured compound; a pale cell with a dashed outline is one nobody has made yet. Click either ' +
      'to inspect it.',
  },
  {
    anchor: () => document.querySelector('.grok-prop-panel'),
    position: ui.hints.POSITION.LEFT,
    title: 'A virtual analog',
    text: 'Select a dashed cell and the Context Panel explains it: the core and substituent were ' +
      'never combined, so its potency is predicted from how each performs elsewhere in the matrix, ' +
      'with the compounds behind every term listed as reference points. Add to make-list collects ' +
      'the ones worth making into the Make list tab.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-control-bar .chem-sar-cart-icon'),
    position: ui.hints.POSITION.LEFT,
    title: 'Collect what is worth making',
    text: 'The cart adds whichever cell is selected in the matrix to the make list, so an analog can ' +
      'be collected without leaving the grid. It takes measured compounds as readily as predicted ' +
      'ones. To collect in bulk, right-click the matrix and add a whole series of predictions at once.',
  },
  {
    anchor: () => Array.from(document.querySelectorAll('.d4-tab-header'))
      .find((e) => e.textContent?.trim() === 'SAR Transfer') ?? null,
    position: ui.hints.POSITION.BOTTOM,
    title: 'Carry the SAR across scaffolds',
    text: 'The SAR Transfer tab pairs cores whose potency trends run in parallel over the substituents ' +
      'they share. Where the correlation holds, an optimization found on one scaffold is expected to ' +
      'carry to the other, and the analogs it argues for are marked in the matrix.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-xfer-panel .chem-sar-grid-host') ??
      Array.from(document.querySelectorAll('.d4-tab-header'))
        .find((e) => e.textContent?.trim() === 'SAR Transfer') ?? null,
    position: ui.hints.POSITION.LEFT,
    title: 'Analogs the transfer argues for',
    text: 'The two series are laid out over the substituents they share, and the columns captioned ' +
      '"predicted" are the ones only one of them has tried. Those cells are dashed, exactly as in the ' +
      'matrix: click one and the Context Panel gives the same Free-Wilson breakdown and the same Add ' +
      'to make-list, so an analog suggested by a transfer joins the same synthesis list.',
  },
  {
    anchor: () => Array.from(document.querySelectorAll('.d4-tab-header'))
      .find((e) => e.textContent?.trim() === 'Make list') ?? null,
    position: ui.hints.POSITION.BOTTOM,
    title: 'The make-list',
    text: 'Everything collected lands here: each analog with its structure, whether it is predicted or ' +
      'already made, its potency and the activity that potency is read on, how much evidence stood ' +
      'behind it, and the series, core and substituent it came from. Open as table hands out a copy ' +
      'to save, export or join; Clear empties the list.',
  },
];

function showSarHint(i: number, previous?: HTMLElement): void {
  previous?.remove();
  if (i >= SAR_HINTS.length)
    return;
  const hint = SAR_HINTS[i];
  const anchor = hint.anchor();
  if (!(anchor instanceof HTMLElement)) {
    showSarHint(i + 1);
    return;
  }
  const content = ui.div();
  content.append(ui.h3(hint.title));
  content.append(ui.divText(hint.text));
  const buttonHost = ui.divH([], {style: {justifyContent: 'flex-end'}});
  content.append(buttonHost);
  const popup = ui.hints.addHint(anchor, content, hint.position);
  buttonHost.append(i === SAR_HINTS.length - 1 ? ui.button('Close', () => popup.remove()) :
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
    let mmpViewer: MatchedMolecularPairsViewer | null = null;
    try {
      await awaitCheck(() => {
        for (const v of tv.viewers) {
          console.log(v.type);
          if (v.type === 'Matched Molecular Pairs Analysis') {
            mmpViewer = v as MatchedMolecularPairsViewer;
            return true;
          }
        }
        return false;
      }, '', 20000);
    } catch (e) {};
    mmpViewer!.helpUrl = 'https://raw.githubusercontent.com/datagrok-ai/public/refs/heads/master/help/datagrok/solutions/domains/chem/chem.md#matched-molecular-pairs';
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
