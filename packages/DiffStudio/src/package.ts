/* eslint-disable max-len */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {solveDefault, solveIVP} from './solver-tools';
import {DiffStudio} from './app';
import {getIVP, IVP, getScriptLines, getScriptParams, ARG_INPUT_KEYS, SCRIPTING, STAGE_COL_NAME} from './scripting-tools';

import {getBallFlightSim} from './demo/ball-flight';
import {CONTROL_EXPR, DF_NAME, LOOP} from './constants';
import {UI_TIME} from './ui-constants';

import {ODEs, SolverOptions} from 'diff-grok';
import {Model} from './model';

import utc from 'dayjs/plugin/utc';
import dayjs from 'dayjs';
import {getLookupsInfo} from './utils';

export const _package = new DG.Package();

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

  @grok.decorators.app({
    name: 'Diff Studio',
    description: 'Solver of ordinary differential equations systems',
    browsePath: 'Compute',
  })
  static async runDiffStudio(): Promise<DG.ViewBase> {
    const path = grok.shell.startUri;
    const toSetStartingPath = (path === window.location.href);

    const proxiView = DG.View.create();

    setTimeout(async () => {
      proxiView.close();

      const solver = new DiffStudio(false);

      if (toSetStartingPath)
        solver.setStartingPath(path);

      const view = await solver.runSolverApp();

      if (view !== null)
        grok.shell.addView(view);
    }, UI_TIME.APP_RUN_SOLVING);

    return proxiView;
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
    let path: string;

    if (!DiffStudio.isStartingUriProcessed) {
      DiffStudio.isStartingUriProcessed = true;
      path = grok.shell.startUri;
    } else
      path = window.location.href;

    const proxiView = DG.View.create();

    setTimeout(async () => {
      proxiView.close();
      const solver = new DiffStudio(false, true, true);
      grok.shell.addView(await solver.getFilePreview(file, path));
    }, UI_TIME.PREVIEW_RUN_SOLVING);

    return proxiView;
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

  @grok.decorators.func({
    description: 'Run model with Diff Studio UI',
  })
  static async runModel(model: string,
    @grok.decorators.param({type: 'int'}) inputsTabDockRatio: number,
    @grok.decorators.param({type: 'int'}) graphsDockRatio: number): Promise<void> {
    const diffStudioModel = new Model(model, {
      inputsTabDockRatio: inputsTabDockRatio,
      graphsDockRatio: graphsDockRatio,
    }, '');

    await diffStudioModel.run();
  }

  @grok.decorators.func({
    'meta': {
      'scriptHandler.language': 'ivp',
      'scriptHandler.extensions': 'ivp',
      'scriptHandler.commentStart': '#',
      'scriptHandler.codeEditorMode': 'python',
      'scriptHandler.parserFunction': 'DiffStudio:ivpLanguageParser',
      'scriptHandler.templateScript': '#name: Template\\n#language: ivp\\n#equations:\\n  dy/dt = -y + sin(t) / t\\n\\n#argument: t\\n  initial = 0.01 {min: 0.01; max: 10}\\n  final = 15 {min: 15; max: 150}\\n  step = 0.01 {min: 0.001; max: 0.1}\\n\\n#inits:\\n  y = 0 {min: 0; max: 9}\\n',
      'icon': 'files/icons/package.png'
    },
    'tags': [
      'scriptHandler'
    ]
  })
  static async ivpLanguageHandler(ivpCall: DG.FuncCall): Promise<void> {
    const params = {...ivpCall.inputs};

    const ivp = getIVP((ivpCall.func as DG.Script).script);
    const scriptText = getScriptLines(ivp).filter(line => line.search(/^[\s]*\/\/language:[\s]*ivp/g) < 0).join('\n');
    const jsScript = DG.Script.create(scriptText);

    const jsCall = jsScript.prepare(params);
    await jsCall.call();
    ivpCall.outputs[DF_NAME] = jsCall.outputs[DF_NAME];
  }

  @grok.decorators.func()
  static ivpLanguageParser(@grok.decorators.param({ type: 'string' }) code: string): DG.Script  {
    const ivp = getIVP(code);

    const lookupsOptions = ivp.inputsLookup ? [getLookupsInfo(ivp.inputsLookup) as any] : [];
    const argOptions = ARG_INPUT_KEYS.map((key) => getOptions(key, ivp.arg[key], CONTROL_EXPR.ARG));
    const initsOptions =  [...ivp.inits.entries()].map(([key, val]) => getOptions(key, val, CONTROL_EXPR.INITS));
    const paramsOptions = ivp.params ? [...ivp.params.entries()].map(([key, val]) => getOptions(key, val, CONTROL_EXPR.PARAMS)) : [];
    const loopOptions = ivp.loop ? [getOptions(LOOP.COUNT_NAME, ivp.loop.count, CONTROL_EXPR.LOOP)] : [];

    const inputs: DG.Property[] = [
      ...lookupsOptions, ...argOptions, ...initsOptions, ...paramsOptions, ...loopOptions
    ].map((propOpts) =>{
      const prop = DG.Property.fromOptions(propOpts);
      return prop;
    });

    const dfProp = DG.Property.fromOptions({
      name: DF_NAME,
      type: DG.TYPE.DATA_FRAME,
      caption: ivp.name,
      // @ts-ignore:next-line
      viewer: getViewersSpec(ivp),
    });

    const ivpScript = DG.Script.fromParams(inputs, [dfProp], code);
    return ivpScript;
  }
}

export * from './package.g';
