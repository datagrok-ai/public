import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {initScriptEditor} from './script-editor';
import {IconTool} from './icon-tool';
import {EntityType} from './constants';
import '../css/styles.css';
import {TestManager} from './package-testing';
import {functionSignatureEditor} from './signature-editor/function-signature-editor';
import {addToJSContextCommand, getMinifiedClassNameMap, hasSupportedType, _renderDevPanel} from './dev-panel';
import {_testDetectorsDialog, _testDetectorsStandard} from './utils/test-detectors';
export * from './package.g';
export const _package = new DG.Package();
let minifiedClassNameMap = {};
export let c: DG.FuncCall;

export class PackageFunctions {
  @grok.decorators.func()
  static renderDevPanel(
    @grok.decorators.param({type: 'object'}) ent: EntityType): Promise<DG.Widget> {
    return _renderDevPanel(ent, minifiedClassNameMap);
  }


  @grok.decorators.func({
    meta: {inspectorPanel: 'true'},
    friendlyName: 'DevTools',
  })
  static _makeInspectorPanel(): DG.Widget {
    return DG.Widget.fromRoot(ui.divText('Custom panel from DevTools'));
  }


  @grok.decorators.autostart({description: 'DevTools autostart function'})
  static autostartTools(): void {
    // Dev pane
    minifiedClassNameMap = getMinifiedClassNameMap();

    grok.events.onAccordionConstructed.subscribe((acc: DG.Accordion) => {
      const ent = acc.context;
      if (ent == null || !hasSupportedType(ent, minifiedClassNameMap)) return;
      const devPane = acc.getPane('Dev');
      if (!devPane)
        acc.addPane('Dev', () => ui.wait(async () => (await PackageFunctions.renderDevPanel(ent)).root));
    });

    grok.events.onContextMenu.subscribe((args) => {
      if (args.args.context instanceof DG.Viewer && (args.args.context as DG.Viewer).type !== DG.VIEWER.GRID)
        addToJSContextCommand(args);
    });

    // script editor, signature editor
    grok.events.onViewAdded.subscribe((view) => {
      if (view.type == 'ScriptView' || view.type == 'DataQueryView')
        functionSignatureEditor(view);
      if (view.type == 'ScriptView')
        initScriptEditor(view);
    });
  }


  @grok.decorators.func({description: 'IconTool'})
  static _IconTool(): void {
    grok.shell.newView('Icon Tool', [new IconTool('Icon Tool')]);
  }


  @grok.decorators.app({
    browsePath: 'Admin',
    name: 'Test Manager',
    'top-menu': 'Tools | Dev | Test Manager',
  })
  static async testManager(): Promise<void> {
    c = grok.functions.getCurrentCall();
    const testManager = new TestManager('Test Manager', true);
    await testManager.init();
  }


  @grok.decorators.func({
    name: 'TestDetectors',
    'top-menu': 'Tools | Dev | Test | Detectors...',
  })
  static testDetectors() {
    _testDetectorsDialog();
  }


  @grok.decorators.func({'top-menu': 'Tools | Dev | Test | Detectors Standard'})
  static async TestDetectorsStandard() {
    const detectorsArray = DG.Func.find({meta: {role: DG.FUNC_TYPES.SEM_TYPE_DETECTOR}});
    const df = await _testDetectorsStandard(detectorsArray);
    grok.shell.addTableView(df);
  }

  @grok.decorators.func()
  static async testFunctions(@grok.decorators.param({type: 'map'}) scope: object) : Promise<DG.DataFrame> {
    const functions = DG.Func.find(scope ?? {});
    const testRuns = {};
    const testTime = [];
    for (const f of functions)
      testRuns[f.name] = await testFunc(f);
    const functionColumn = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'function', Object.keys(testRuns));
    const results = Object.values(testRuns).map((runs) => {
      let totalTime = 0;
      if (!runs) {
        testTime.push(totalTime);
        return null;
      }
      const failedTests = {};
      const cntTotal = Object.keys(runs).length;
      let cntPassed = 0;
      for (const [test, result] of Object.entries(runs)) {
        totalTime += result.time;
        if (result.output == true)
          cntPassed++;
        else
          failedTests[test] = result.output;
      }
      testTime.push(totalTime);
      return cntPassed === cntTotal ? 'OK' : `${cntPassed} / ${cntTotal}. Failed: ${JSON.stringify(failedTests)}`;
    });
    const resultColumn = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'result', results);
    const timeColumn = DG.Column.fromList(DG.COLUMN_TYPE.INT, 'ms', testTime);
    return DG.DataFrame.fromColumns([functionColumn, resultColumn, timeColumn]);
  }

  @grok.decorators.func({meta: {vectorFunc: 'true'}})
  static async testFunction(
    @grok.decorators.param({type: 'column<int>'}) col1: DG.Column,
    @grok.decorators.param({type: 'column<string>'}) col2: DG.Column,
    @grok.decorators.param({type: 'column<double>'}) col3: DG.Column,
    @grok.decorators.param({type: 'list<string>', options: {optional: true}}) out?: string[]): Promise<DG.DataFrame> {
    function createNewCol1() {
      const res1 = DG.Column.int('newCol1', col1.length);
      res1.init((i) => col1.getNumber(i) + 1);
      return res1;
    }
    function createNewCol2() {
      const res2 = DG.Column.string('newCol2', col2.length);
      res2.init((i) => col2.get(i) + ' and 123');
      return res2;
    }
    function createNewCol3() {
      const res3 = DG.Column.float('newCol3', col3.length);
      res3.init((i) => col3.getNumber(i) + 5.5);
      return res3;
    }
    const colCreationFuncs: {[colName: string]: () => DG.Column} = {
      'newCol1': createNewCol1,
      'newCol2': createNewCol2,
      'newCol3': createNewCol3,
    };

    console.log(out);
    const colList: DG.Column[] = [];
    if (out == undefined || out.length === 0)
      return DG.DataFrame.fromColumns([createNewCol1(), createNewCol2(), createNewCol3()]);
    else {
      for (const colName of out) {
        if (colCreationFuncs[colName] != undefined)
          colList.push(colCreationFuncs[colName]());
      }
    }
    return colList.length > 0 ? DG.DataFrame.fromColumns(colList) : DG.DataFrame.create(col1.length);
  }

  @grok.decorators.func({
    meta: {vectorFunc: 'true'},
    outputs: [{name: 'result', type: 'dataframe', options: {action: 'join(data)'}}],
  })
  static async testFunctionJoin(
    data: DG.DataFrame,
    @grok.decorators.param({type: 'column<int>'}) col1: DG.Column,
    @grok.decorators.param({type: 'column<string>'}) col2: DG.Column,
    @grok.decorators.param({type: 'column<double>'}) col3: DG.Column,
    @grok.decorators.param({type: 'list<string>', options: {optional: true}}) out?: string[]): Promise<DG.DataFrame> {
    function createNewCol1() {
      const res1 = DG.Column.int('joinedCol1', col1.length);
      res1.init((i) => col1.getNumber(i) + 1);
      return res1;
    }
    function createNewCol2() {
      const res2 = DG.Column.string('joinedCol2', col2.length);
      res2.init((i) => col2.get(i) + ' joined');
      return res2;
    }
    function createNewCol3() {
      const res3 = DG.Column.float('joinedCol3', col3.length);
      res3.init((i) => col3.getNumber(i) + 10.5);
      return res3;
    }
    const colCreationFuncs: {[colName: string]: () => DG.Column} = {
      'joinedCol1': createNewCol1,
      'joinedCol2': createNewCol2,
      'joinedCol3': createNewCol3,
    };

    const colList: DG.Column[] = [];
    if (out == undefined || out.length === 0)
      return DG.DataFrame.fromColumns([createNewCol1(), createNewCol2(), createNewCol3()]);
    else {
      for (const colName of out) {
        if (colCreationFuncs[colName] != undefined)
          colList.push(colCreationFuncs[colName]());
      }
    }
    return colList.length > 0 ? DG.DataFrame.fromColumns(colList) : DG.DataFrame.create(col1.length);
  }


  @grok.decorators.func({
    name: 'ExceptionFunc',
    outputs: [{name: 'result', type: 'int'}],
  })
  static exceptionFunc(
    @grok.decorators.param({type: 'int'}) a: number): number {
    if (a === 0)
      throw 'exception';
    else
      a++;
    return a;
  }
}


export async function testFunc(f: DG.Func): Promise<{[key: string]: {output: any, time: number}}> {
  const tests = f.options['test'];
  if (tests == null || Array.isArray(tests) && !tests.length)
    return null;
  const results = {};
  function addNamespace(s: string): string {
    return s.replace(new RegExp(f.name, 'gi'), f.nqName);
  }
  for (const test of tests) {
    const start = Date.now();
    const output = await grok.functions.eval(addNamespace(test));
    const time = Date.now() - start;
    results[test] = {output, time};
  }
  return results;
}
