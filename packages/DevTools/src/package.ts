import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {initScriptEditor} from './script-editor';
import {IconTool} from './icon-tool';
import {EntityType} from './constants';
import '../css/styles.css';
// import * as tests from './tests/test-examples';
import {TestManager} from './package-testing';
import {functionSignatureEditor} from './signature-editor/function-signature-editor';
import {addToJSContextCommand, getMinifiedClassNameMap, hasSupportedType, _renderDevPanel} from './dev-panel';
import {_testDetectorsDialog, _testDetectorsStandard} from './utils/test-detectors';

export const _package = new DG.Package();
let minifiedClassNameMap = {};
export let c: DG.FuncCall;

//name: renderDevPanel
//tags: dev-tools
//input: object ent
//output: widget panel
export function renderDevPanel(ent: EntityType): Promise<DG.Widget> {
  return _renderDevPanel(ent, minifiedClassNameMap);
}

//friendlyName: DevTools
//output: widget panel
//meta.inspectorPanel: true
export function _makeInspectorPanel(): DG.Widget {
  return DG.Widget.fromRoot(ui.divText('Custom panel from DevTools'));
}

//description: DevTools autostart function
//tags: autostart
export function autostartTools(): void {
  // Dev pane
  minifiedClassNameMap = getMinifiedClassNameMap();

  grok.events.onAccordionConstructed.subscribe((acc: DG.Accordion) => {
    const ent = acc.context;
    if (ent == null || !hasSupportedType(ent, minifiedClassNameMap)) return;
    const devPane = acc.getPane('Dev');
    if (!devPane)
      acc.addPane('Dev', () => ui.wait(async () => (await renderDevPanel(ent)).root));
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

//description: IconTool
export function _IconTool(): void {
  grok.shell.newView('Icon Tool', [new IconTool('Icon Tool')]);
}

//name: Test Manager
//top-menu: Tools | Dev | Test Manager
//meta.browsePath: Admin
//tags: app
export async function testManager(): Promise<void> {
  c = grok.functions.getCurrentCall();
  const testManager = new TestManager('Test Manager', true);
  await testManager.init();
}

//name: TestDetectors
//top-menu: Tools | Dev | Test | Detectors...
export function testDetectors() {
  _testDetectorsDialog();
}

//name: TestDetectorsStandard
//top-menu: Tools | Dev | Test | Detectors Standard
export async function TestDetectorsStandard() {
  const detectorsArray = DG.Func.find({tags: ['semTypeDetector']});
  const df = await _testDetectorsStandard(detectorsArray);
  grok.shell.addTableView(df);
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

//name: testFunctions
//input: map scope
//output: dataframe result
export async function testFunctions(scope: object) {
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

//name: testFunction
//description: test func
//meta.vectorFunc: true
//input: column<int> col1
//input: column<string> col2
//input: column<double> col3
//input: list<string> out { optional: true }
//output: dataframe result
export async function testFunction(col1: DG.Column, col2: DG.Column, col3: DG.Column, out?: string[]): Promise<DG.DataFrame> {
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

//name: ExceptionFunc
//language: javascript
//input: int a
//output: int out
export function exceptionFunc(a: number): number {
  if (a === 0)
    throw 'exception';
  else
    a++;
  return a;
}
