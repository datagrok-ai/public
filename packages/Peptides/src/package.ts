/* eslint-disable max-len */
/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {View} from 'datagrok-api/dg';

import {analyzePeptidesUI} from './widgets/peptides';
import {manualAlignmentWidget} from './widgets/manual-alignment';
import {MonomerPosition, MostPotentResidues} from './viewers/sar-viewer';
import {getTreeHelper, ITreeHelper} from '@datagrok-libraries/bio/src/trees/tree-helper';
import {LogoSummaryTable} from './viewers/logo-summary';
import {MonomerWorks} from '@datagrok-libraries/bio/src/monomer-works/monomer-works';
import {PeptidesModel} from './model';
import {macromoleculeSarFastaDemoUI} from './demo/fasta';
import {u2} from '@datagrok-libraries/utils/src/u2';
import {ClusterMaxActivityViewer} from './viewers/cluster-max-activity-viewer';
import {LSTPieChartRenderer} from './utils/cell-renderer';
import {PeptideUtils} from './peptideUtils';
import {SequencePositionStatsViewer} from './viewers/position-statistics-viewer';
import {MutationCliffsViewer} from './viewers/mutation-cliffs-viewer';

let monomerWorks: MonomerWorks | null = null;
let treeHelper: ITreeHelper;
export const _package = new DG.Package();
export * from './package.g';

/** Temporary polyfill */

function getDecoratorFunc() {
  return function(args: any) {
    return function(
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor,
    ) { };
  };
}

// Ensure decorators object exists and polyfill missing decorators
if (!grok.decorators)
  (grok as any).decorators = {};


const decorators = [
  'func', 'init', 'param', 'panel', 'editor', 'demo', 'app',
  'appTreeBrowser', 'fileHandler', 'fileExporter', 'model', 'viewer', 'filter', 'cellRenderer', 'autostart',
  'dashboard', 'folderViewer', 'semTypeDetector', 'packageSettingsEditor', 'functionAnalysis', 'converter',
  'fileViewer', 'model', 'treeBrowser', 'polyfill',
];

decorators.forEach((decorator) => {
  if (!(grok.decorators as any)[decorator])
    (grok.decorators as any)[decorator] = getDecoratorFunc();
});

/** End temporary polyfill */

export function getMonomerWorksInstance(): MonomerWorks | null {
  return monomerWorks;
}


export function getTreeHelperInstance(): ITreeHelper {
  return treeHelper;
}


async function openDemoData(chosenFile: string): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('Loading Peptides');
  const path = _package.webRoot + 'files/' + chosenFile;
  const peptides = await grok.data.loadTable(path);
  peptides.name = 'Peptides';
  const view = grok.shell.addTableView(peptides);
  view.name = 'PeptidesView';
  grok.shell.windows.showProperties = true;
  pi.close();
}


export class PackageFunctions {
  @grok.decorators.init({tags: ['init']})
  static async initPeptides(): Promise<void> {
    try {
      monomerWorks ??= new MonomerWorks(await grok.functions.call('Bio:getBioLib'));
      treeHelper ??= await getTreeHelper();
      await PeptideUtils.loadComponents();
    } catch (e) {
      grok.log.error(e as string);
    }
  }

  @grok.decorators.func()
  static Peptides(): DG.View {
    const appHeader = u2.appHeader({
      iconPath: _package.getIconUrl(),
      learnMoreUrl: 'https://github.com/datagrok-ai/public/blob/master/help/domains/bio/peptides.md',
      description:
        '- Automatically recognizes peptides in your data\n' +
        '- Invariant map and mutation cliffs\n' +
        '- Logo plots to explore sequence composition\n' +
        '- Hierarchical clustering\n' +
        '- Sequence space to analyze clustering and activity cliffs\n' +
        '- Finds statistically significant changes in activity for monomer/positions\n',
    });

    const windows = grok.shell.windows;
    windows.showToolbox = false;
    windows.showHelp = false;
    windows.showProperties = false;

    const view = View.create();
    view.name = 'Peptides';
    ui.appendAll(view.root, [
      appHeader,
      ui.divH([
        ui.button('Simple demo', () => openDemoData('aligned.csv'), ''),
        ui.button('Complex demo', () => openDemoData('aligned_2.csv'), ''),
        ui.button('HELM demo', () => openDemoData('aligned_3.csv'), ''),
      ]),
    ]);
    return view;
  }


  @grok.decorators.func({
    'top-menu': 'Bio | Analyze | SAR...',
    'name': 'Bio Peptides',
    'outputs': [],
  })
  static peptidesDialog(): DG.Dialog | null {
    if (!grok.shell.t || !grok.shell.t.columns.bySemType('Macromolecule')?.length) {
      grok.shell.warning('SAR Analysis requires an active table with Macromolecule column');
      return null;
    }

    if (!DG.Utils.firstOrNull(grok.shell.t.columns.numerical)) {
      grok.shell.warning('SAR Analysis requires an active table with at least one numerical column for activity');
      return null;
    }

    const analyzeObject = analyzePeptidesUI(grok.shell.t);
    const dialog = ui.dialog('Analyze Peptides').add(analyzeObject.host).onOK(async () => {
      const startSuccess = analyzeObject.callback();
      if (!startSuccess)
        dialog.show();
    });
    return dialog.show();
  }


  @grok.decorators.func()
  static async testInitFunctionPeptides(
    @grok.decorators.param({'type': 'viewer'}) v: DG.Viewer): Promise<void> {
    grok.shell.info('Test init function for Peptides package');
    grok.shell.info('Viewer name: ' + v.dataFrame.name);
    await new Promise<void>((r) => setTimeout(r, 1000));
  }


  @grok.decorators.panel({
    name: 'Peptides',
    meta: {role: 'widgets'},
    tags: ['widgets', 'panel'],
  })
  static peptidesPanel(
    @grok.decorators.param({'options': {'semType': 'Macromolecule'}}) col: DG.Column): DG.Widget {
    if (!col.dataFrame || !DG.Utils.firstOrNull(col.dataFrame.columns.numerical))
      return new DG.Widget(ui.divText('SAR Analysis requires an active table with at least one numerical column for activity'));
    const analyzeObject = analyzePeptidesUI(col.dataFrame, col);
    return new DG.Widget(analyzeObject.host);
  }


  @grok.decorators.func({
    meta: {icon: 'files/icons/peptide-sar-viewer.svg', role: 'viewer'},
    name: 'Sequence Variability Map',
    tags: ['viewer'],
    description: 'Peptides Sequence Variability Map Viewer',
    outputs: [{type: 'viewer', name: 'result'}],
  })
  static monomerPosition(): DG.Viewer {
    return new MonomerPosition();
  }


  @grok.decorators.func({
    meta: {icon: 'files/icons/peptide-sar-vertical-viewer.svg', role: 'viewer'},
    name: 'Most Potent Residues',
    tags: ['viewer'],
    description: 'Peptides Most Potent Residues Viewer',
    outputs: [{type: 'viewer', name: 'result'}],
  })
  static mostPotentResidues(): DG.Viewer {
    return new MostPotentResidues();
  }

  @grok.decorators.func({
    meta: {icon: 'files/icons/sequence-statistics-viewer.svg', role: 'viewer'},
    name: 'Sequence Mutation Cliffs',
    description: 'Mutation Cliffs Line Chart',
    tags: ['viewer'],
    outputs: [{type: 'viewer', name: 'result'}],
  })
  static mutationCliffs(): DG.Viewer {
    return new MutationCliffsViewer();
  }


  @grok.decorators.func({
    meta: {icon: 'files/icons/logo-summary-viewer.svg', role: 'viewer'},
    name: 'Logo Summary Table',
    tags: ['viewer'],
    outputs: [{type: 'viewer', name: 'result'}],
  })
  static logoSummaryTable(): DG.Viewer {
    return new LogoSummaryTable();
  }


  @grok.decorators.func({
    meta: {icon: 'files/icons/sequence-statistics-viewer.svg', role: 'viewer'},
    name: 'Sequence Position Statistics',
    outputs: [{type: 'viewer', name: 'result'}],
    tags: ['viewer'],
  })
  static sequencePositionStatistics(): DG.Viewer {
    return new SequencePositionStatsViewer();
  }


  @grok.decorators.func({
    name: 'Active peptide selection',
    outputs: [{type: 'viewer', name: 'result'}],
    meta: {role: 'viewer'},
    tags: ['viewer'],
  })
  static clusterMaxActivity(): DG.Viewer {
    return new ClusterMaxActivityViewer();
  }


  @grok.decorators.panel({
    name: 'Manual Alignment',
    meta: {role: 'widgets'},
    tags: ['widgets', 'panel'],
  })
  static manualAlignment(
    @grok.decorators.param({'options': {'semType': 'Monomer'}}) _monomer: string): DG.Widget {
    //TODO: recalculate Molfile and Molecule panels on sequence update
    const df = grok.shell.t;
    const model: PeptidesModel | null = df?.temp[PeptidesModel.modelName];
    if (!model)
      return new DG.Widget(ui.divText('Manual alignment works with peptides analysis'));


    const col = df.getCol(model.settings!.sequenceColumnName!);
    return manualAlignmentWidget(col, df);
  }


  @grok.decorators.func({
    meta: {
      demoPath: 'Bioinformatics | Peptide SAR',
      isDemoDashboard: 'true',
    },
    name: 'Peptide SAR',
    description: 'Peptide SAR Analysis demo on peptide sequences in FASTA format',
  })
  static async macromoleculeSarFastaDemo(): Promise<void> {
    return await macromoleculeSarFastaDemoUI();
  }


  @grok.decorators.func({
    meta: {
      cellType: 'lst-pie-chart',
      gridChart: 'true',
      role: 'cellRenderer',
    },
    tags: ['cellRenderer'],
    name: 'LST Pie Chart',
    outputs: [{type: 'grid_cell_renderer', name: 'result'}],
  })
  static lstPiechartCellRenderer(): DG.GridCellRenderer {
    return new LSTPieChartRenderer();
  }
}
