import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {scriptEditor} from './script-editor';
import {IconTool} from './icon-tool';
import {EntityType} from './constants';
import '../css/styles.css';
import * as tests from './tests/test-examples';
import {TestManager} from './package-testing';
import {functionSignatureEditor} from './function-signature-editor';
import {addToJSContextCommand, getMinifiedClassNameMap, _renderDevPanel} from './dev-panel';
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


//tags: autostart
export function describeCurrentObj(): void {
  minifiedClassNameMap = getMinifiedClassNameMap();

  grok.events.onAccordionConstructed.subscribe((acc: DG.Accordion) => {
    const ent = acc.context;
    if (ent == null) return;
    const devPane = acc.getPane('Dev');
    if (!devPane)
      acc.addPane('Dev', () => ui.wait(async () => (await renderDevPanel(ent)).root));
  });

  grok.events.onContextMenu.subscribe((args) => {
    if (args.args.context instanceof DG.Viewer)
      addToJSContextCommand(args);
  });
}

//description: ScriptEditor
//tags: autostart
export function _scriptEditor(): void {
  grok.events.onViewAdded.subscribe((view) => {
    if (view.type == 'ScriptView')
      scriptEditor(view);
  });
}

//description: FunctionSignatureEditor
//tags: autostart
export function _functionSignatureEditor(): void {
  grok.events.onViewAdded.subscribe((view) => {
    if (view.type == 'ScriptView')
      functionSignatureEditor(view);
  });
}

//description: IconTool
export function _IconTool(): void {
  grok.shell.newView('Icon Tool', [new IconTool('Icon Tool')]);
}


//name: TestManager
//top-menu: Tools | Dev | Test Manager
//tags: app
export async function testManager(): Promise<void> {
  c = grok.functions.getCurrentCall();
  const testManager = new TestManager('Test Manager');
  await testManager.init();
}


//name: TestDetectors
//top-menu: Tools | Dev | Test Detectors
export function testDetectors() {
  _testDetectorsDialog();
}

//name: TestDetectorsStandard
//top-menu: Tools | Dev | Test Detectors Standard
export async function TestDetectorsStandard() {
  const detectorsArray = DG.Func.find({tags: ['semTypeDetector']});
  const df = await _testDetectorsStandard(detectorsArray);
  grok.shell.addTableView(df);
}

async function testFunc(f: DG.Func): Promise<{[key: string]: any}> {
  const tests = f.options['test'];
  if (tests == null || Array.isArray(tests) && !tests.length)
    return null;

  const results = {};

  function addNamespace(s: string): string {
    return s.replace(new RegExp(f.name, 'gi'), f.nqName);
  }

  for (const test of tests) {
    results[test] = await grok.functions.eval(addNamespace(test));
  }

  return results;
}

//name: testFunctions
//input: string scope = "" [JSON string with the search filter for functions]
//output: dataframe result
export async function testFunctions(scope: string = '') {
  const functions = DG.Func.find(scope ? JSON.parse(scope) : {package: _package.name});
  const testRuns = {};
  for (const f of functions) {
    testRuns[f.name] = await testFunc(f);
  }
  const functionColumn = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'function', Object.keys(testRuns));
  const results = Object.values(testRuns).map((runs) => {
    if (!runs)
      return null;
    const failedTests = {};
    const cntTotal = Object.keys(runs).length;
    let cntPassed = 0;
    for (const [test, result] of Object.entries(runs)) {
      if (result == true)
        cntPassed++;
      else
        failedTests[test] = result;
    }
    return cntPassed === cntTotal ? 'OK' : `${cntPassed} / ${cntTotal}. Failed: ${JSON.stringify(failedTests)}`;
  });
  const resultColumn = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'result', results);
  return DG.DataFrame.fromColumns([functionColumn, resultColumn]);
}

