/* eslint-disable max-len */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {solveDefault, solveIVP} from './solver-tools';
import {DiffStudio} from './app';
import {getIVP, IVP, getScriptLines, getScriptParams, ARG_INPUT_KEYS, SCRIPTING} from './scripting-tools';

import {getBallFlightSim} from './demo/ball-flight';
import {PK_PD_DEMO} from './demo/pk-pd';
import {BIOREACTOR_DEMO} from './demo/bioreactor';
import {CONTROL_EXPR, DF_NAME} from './constants';
import {UI_TIME} from './ui-constants';

import {ODEs, SolverOptions} from '@datagrok/diff-grok';
import {Model} from './model';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

export class PackageFunctions {
  @grok.decorators.init({})
  static async init() {}

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
        options: {caption: 'Trajectory', viewer: 'Line chart(block: 60, multiAxis: "false", multiAxisLegendPosition: "RightCenter", autoLayout: "false", showAggrSelectors: "false") | Grid(block: 40)'},
      },
    ],
  })
  static ballFlight(
    @grok.decorators.param({options: {initialValue: '0.01', category: 'Ball', caption: 'Diameter', units: 'm', min: '0.01', max: '0.3'}}) dB: number,
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

  @grok.decorators.model({
    name: 'PK-PD',
    description: 'In-browser two-compartment pharmacokinetic-pharmacodynamic (PK-PD) simulation',
    icon: 'files/icons/pkpd.png',
  })
  static async pkPdNew(): Promise<void> {
    await PK_PD_DEMO.run();
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
    await PK_PD_DEMO.runDemo();
  }

  @grok.decorators.model({
    description: 'Controlled fab-arm exchange mechanism simulation',
    icon: 'files/icons/bioreactor.png',
  })
  static async Bioreactor(): Promise<void> {
    await BIOREACTOR_DEMO.run();
  }

  @grok.decorators.demo({
    name: 'Bioreactor Demo',
    description: 'In-browser simulation of controlled fab-arm exchange mechanism',
    demoPath: 'Compute | Bioreactor',
    test: {test: 'demoBioreactor()', wait: '100'},
  })
  static async demoBioreactor(): Promise<any> {
    await BIOREACTOR_DEMO.runDemo();
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

    const argOptions = ARG_INPUT_KEYS.map((key) => DiffStudio.getOptions(key, ivp.arg[key], CONTROL_EXPR.ARG));
    const initsOptions =  [...ivp.inits.entries()].map(([key, val]) => DiffStudio.getOptions(key, val, CONTROL_EXPR.INITS));
    const paramsOptions = ivp.params ? [...ivp.params.entries()].map(([key, val]) => DiffStudio.getOptions(key, val, CONTROL_EXPR.PARAMS)) : [];
    const loopOptions = ivp.loop ? [DiffStudio.getOptions(SCRIPTING.COUNT, ivp.loop.count, CONTROL_EXPR.LOOP)] : [];

    const inputs: DG.Property[] = [
      ...argOptions, ...initsOptions, ...paramsOptions, ...loopOptions
    ].map((propOpts) =>{
      const prop = DG.Property.fromOptions(propOpts);
      return prop;
    });

    const dfProp = DG.Property.fromOptions({
      name: DF_NAME,
      type: DG.TYPE.DATA_FRAME,
    });

    dfProp.options.viewer = 'Grid() | Line chart()';

    const ivpScript = DG.Script.fromParams(inputs, [dfProp], code);
    return ivpScript;
  }
}

export * from './package.g';
