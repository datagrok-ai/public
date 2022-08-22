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

