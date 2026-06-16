/* eslint-disable max-len */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {solveDefault, solveIVP} from './solver-tools';
import {DiffStudio} from './app';
import {DiffStudioFacetViewer} from './diff-studio-facet-viewer';
import {DiffStudioHub} from './hub';
import {getIVP, IVP, getScriptLines, getScriptParams} from './scripting-tools';

import {getBallFlightSim} from './demo/ball-flight';

export {Model} from './model';

import {DF_NAME} from './constants';
import {PATH, TITLE, UI_TIME} from './ui-constants';

import {ODEs, SolverOptions} from 'diff-grok';
import {Model, ModelInfo} from './model';

import utc from 'dayjs/plugin/utc';
import dayjs from 'dayjs';

export const _package = new DG.Package();

/** Default demo layout for `.ivp`-driven demos. */
const DEMO_UI_OPTS = {inputsTabDockRatio: 0.17, graphsDockRatio: 0.85};

/** Run a model shipped under `files/<modelFile>` as a Diff Studio demo — equations come from the
 *  `.ivp` file (single source of truth). Help text is the companion `<helpFile>` readme when given,
 *  otherwise it falls back to text derived from the model's `#description`. */
async function runIvpDemo(modelFile: string, helpFile?: string): Promise<void> {
  const equations = await _package.files.readAsText(modelFile);
  let info = '';
  if (helpFile) {
    try {
      info = await _package.files.readAsText(helpFile);
    } catch {
      // readme missing — fall back to #description below
    }
  }
  if (!info) {
    const descr = (/^#description:\s*(.+)$/m.exec(equations)?.[1] ?? '').trim();
    info = `# Model\n${descr}\n\n# Try\nInteractive results update as you change the inputs.`;
  }
  await new Model({equations, uiOptions: DEMO_UI_OPTS, info}).runDemo();
}

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

export class PackageFunctions {
  @grok.decorators.init({})
  static async init() {
    dayjs.extend(utc);
  }

  @grok.decorators.func({})
  static dock(): void {
    const df = grok.data.demo.demog(100);
    grok.shell.windows.showToolbox = false;
    grok.shell.windows.showBrowse = true;
    const view = grok.shell.addTableView(df);

    setTimeout(() => {
      const chart = DG.Viewer.boxPlot(df);
      let node = view.dockManager.dock(chart, 'right', null, 'ANOVA');

      const div = ui.div([
        ui.button('Click me', () => {
          grok.shell.info('Clicked');
        }),
      ]);
      node = view.dockManager.dock(div, 'down', node, 'Title', 0.3);

      const grid = DG.Viewer.grid(grok.data.demo.demog(100));
      view.dockManager.dock(grid, 'fill', node);
    }, 100);
  }

  @grok.decorators.func({})
  static solve(@grok.decorators.param({type: 'object'}) problem: ODEs): DG.DataFrame {
    return solveDefault(problem);
  }

  @grok.decorators.func({})
  static solveEquations(
    @grok.decorators.param({type: 'object'}) problem: ODEs,
    @grok.decorators.param({type: 'object'}) options: Partial<SolverOptions>): DG.DataFrame {
    return solveIVP(problem, options);
  }

  @grok.decorators.func({
    name: 'DiffStudio Facet',
    description: 'Faceted grid of line charts, one per output variable, for Diff Studio solutions',
    meta: {showInGallery: 'false', role: 'viewer'},
    outputs: [{type: 'viewer', name: 'result'}],
  })
  static diffStudioFacetViewer(): DG.JsViewer {
    return new DiffStudioFacetViewer();
  }

  @grok.decorators.app({
    name: 'Diff Studio',
    description: 'Solver of ordinary differential equations systems',
    browsePath: 'Compute',
  })
  static async runDiffStudio(): Promise<DG.ViewBase> {
    const path = grok.shell.startUri;

    const wasProcessed = DiffStudio.isStartingUriProcessed;
    DiffStudio.isStartingUriProcessed = true;

    const isDeepLink = path.includes(PATH.MODEL) ||
      path.includes(PATH.PARAM) ||
      path.includes(`/${TITLE.TEMPL}`) ||
      path.includes(`/${TITLE.LIBRARY}`) ||
      path.includes(`/${TITLE.RECENT}`);

    if (wasProcessed || !isDeepLink) {
      // Dedup: if a hub view is already open (stamped on `view.temp`),
      // activate it instead of creating a duplicate.
      const existing = DiffStudio.findHubView();
      if (existing) {
        grok.shell.v = existing;
        return existing;
      }
      const hub = new DiffStudioHub();
      DiffStudio.stampHubView(hub.view);
      DiffStudio.currentHubRenderer = () => hub.render();
      const sub = grok.events.onViewRemoved.subscribe((removed: DG.View) => {
        if (removed.temp?.[DiffStudio.HUB_TAG] === true) {
          DiffStudio.currentHubRenderer = null;
          sub.unsubscribe();
        }
      });
      hub.renderHeader();
      setTimeout(() => hub.renderRest(), 0);
      return hub.view;
    }

    const toSetStartingPath = (path === window.location.href);
    // Return an empty proxy view synchronously and mount the real solver view
    // asynchronously via grok.shell.addView. This is required for correct
    // integration with the platform's app-launching flow — DO NOT change this
    // to return the actual solver view directly.
    const proxyView = DG.View.create();

    setTimeout(async () => {
      proxyView.close();

      const solver = new DiffStudio(false);

      if (toSetStartingPath)
        solver.setStartingPath(path);

      const view = await solver.runSolverApp();

      if (view !== null) {
        grok.shell.windows.showToolbox = false;
        grok.shell.windows.showBrowse = true;
        grok.shell.addView(view);
      }
    }, UI_TIME.APP_RUN_SOLVING);

    return proxyView;
  }

  @grok.decorators.demo({
    name: 'Diff Studio Demo',
    description: 'Interactive solver of ordinary differential equations (ODE)',
    test: {test: 'runDiffStudioDemo()', wait: '100'},
    demoPath: 'Compute | Diff Studio',
  })
  static async runDiffStudioDemo(): Promise<void> {
    const solver = new DiffStudio();
    await solver.runSolverDemoApp();
  }

  @grok.decorators.fileHandler({
    ext: 'ipv',
  })
  static async ivpFileHandler(content: string) {
    const solver = new DiffStudio();
    await solver.runSolverApp(content);

    return [];
  }

  @grok.decorators.fileViewer({fileViewer: 'ivp'})
  static async previewIvp(file: DG.FileInfo): Promise<DG.View> {
    const wasProcessed = DiffStudio.isStartingUriProcessed;
    DiffStudio.isStartingUriProcessed = true;

    // Only honor URL params on the initial deep-link load. On subsequent file clicks,
    // window.location.href still holds the previous file's params and would be misapplied.
    const path = wasProcessed ? file.fullPath : grok.shell.startUri;

    // Return an empty proxy view synchronously and mount the real preview
    // asynchronously via grok.shell.addView. This is required for correct
    // integration with the platform's file-viewer flow — DO NOT change this
    // to return the actual preview view directly.
    const proxyView = DG.View.create();

    setTimeout(async () => {
      proxyView.close();
      const solver = new DiffStudio(false, true, true);
      grok.shell.windows.showToolbox = false;
      grok.shell.windows.showBrowse = true;
      grok.shell.addView(await solver.getFilePreview(file, path));
    }, UI_TIME.PREVIEW_RUN_SOLVING);

    return proxyView;
  }

  @grok.decorators.func()
  static async runDiffStudioTreeBrowser(treeNode: DG.TreeViewGroup) {
    new DiffStudio(false, false, false, {treeNode: treeNode, browsePanel: grok.shell.browsePanel});
  }

  @grok.decorators.model({
    name: 'Ball flight',
    description: 'Ball flight simulation',
    editor: 'Compute2:RichFunctionViewEditor',
    runOnOpen: 'true',
    runOnInput: 'true',
    features: '{"sens-analysis": true, "fitting": true}',
    icon: 'files/icons/ball.png',
    // @ts-expect-error
    help: 'ball-flight.md',
    dockSpawnConfig: '{"Trajectory / Grid": {"dock-spawn-dock-ratio": 0.3, "dock-spawn-dock-type": "right", "dock-spawn-dock-to": "Trajectory / Line chart"}, "Output": {"dock-spawn-dock-ratio": 0.15, "dock-spawn-dock-type": "down", "dock-spawn-dock-to": "Trajectory / Line chart"}}',
    outputs: [
      {
        name: 'maxDist',
        type: 'double',
        options: {caption: 'Max distance'},
      },
      {
        name: 'maxHeight',
        type: 'double',
        options: {caption: 'Max height'},
      },
      {
        name: 'df',
        type: 'dataframe',
        options: {caption: 'Trajectory', viewer: 'Line chart(multiAxis: "false", multiAxisLegendPosition: "RightCenter", autoLayout: "false", showAggrSelectors: "false") | Grid()'},
      },
    ],
  })
  static ballFlight(
    // @ts-expect-error
    @grok.decorators.param({options: {initialValue: '0.01', category: 'Ball', caption: 'Diameter', units: 'm', min: '0.01', max: '0.3', minFormula: 'roB / 20000', maxFormula: 'roB / 4000'}}) dB: number,
    @grok.decorators.param({options: {initialValue: '200', category: 'Ball', caption: 'Density', description: 'Material density', units: 'kg/m^3', min: '200', max: '1200'}}) roB: number,
    @grok.decorators.param({options: {initialValue: '50', category: 'Throw parameters', caption: 'Velocity', min: '40', max: '60', units: 'm/sec'}}) v: number,
    @grok.decorators.param({options: {initialValue: '45', category: 'Throw parameters', caption: 'Angle', min: '20', max: '70', units: 'deg'}}) a: number) {
    const simlulation = getBallFlightSim(v, Math.PI * a / 180, dB, roB);
    return {
      df: simlulation,
      maxDist: simlulation.col('Distance').stats.max,
      maxHeight: simlulation.col('Height').stats.max,
    };
  }

  @grok.decorators.func({
    description: 'Return serialized initial value problem for ordinary differential equations',
    outputs: [{name: 'serialization', type: 'object'}],
  })
  static serializeEquations(problem: string): IVP {
    return getIVP(problem);
  }

  @grok.decorators.func({
    description: 'Perform ODEs serialization to JS-code',
  })
  static odesToCode(serialization: IVP): string {
    return getScriptLines(serialization).join('\n');
  }

  @grok.decorators.func({
    description: 'Solve initial value problem for ordinary differential equations',
  })
  static async solveODE(problem: string): Promise<DG.DataFrame> {
    const ivp = getIVP(problem);
    const code = getScriptLines(ivp).join('\n');

    const script = DG.Script.create(code);
    const params = getScriptParams(ivp);
    const call = script.prepare(params);

    await call.call();

    return call.outputs[DF_NAME];
  }

  @grok.decorators.demo({
    name: 'PK-PD Simulation Demo',
    description: 'In-browser two-compartment pharmacokinetic-pharmacodynamic (PK-PD) simulation',
    demoPath: 'Compute | PK-PD Modeling',
    test: {
      test: 'demoSimPKPD()',
      wait: '100',
    },
  })
  static async demoSimPKPD(): Promise<any> {
    await runIvpDemo('models/pk-pd.ivp', 'models/pk-pd.md');
  }


  @grok.decorators.demo({
    name: 'Bioreactor Demo',
    description: 'In-browser simulation of controlled fab-arm exchange mechanism',
    demoPath: 'Compute | Bioreactor',
    test: {test: 'demoBioreactor()', wait: '100'},
  })
  static async demoBioreactor(): Promise<any> {
    await runIvpDemo('models/bioreactor.ivp', 'models/bioreactor.md');
  }


  @grok.decorators.func({
    description: 'Run model with Diff Studio UI',
  })
  static async runModel(model: string,
    @grok.decorators.param({type: 'int'}) inputsTabDockRatio: number,
    @grok.decorators.param({type: 'int'}) graphsDockRatio: number): Promise<void> {
    const modelInfo: ModelInfo = {
      equations: model,
      uiOptions: {
        inputsTabDockRatio: inputsTabDockRatio,
        graphsDockRatio: graphsDockRatio,
      },
      info: '',
    };

    const diffStudioModel = new Model(modelInfo);

    await diffStudioModel.run();
  }
}

export * from './package.g';
