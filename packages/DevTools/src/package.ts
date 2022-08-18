import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {scriptEditor} from './script-editor';
import {IconTool} from './icon-tool';
import {EntityType} from './constants';
import '../css/styles.css';
import * as tests from './tests/test-examples';
import {testManagerView, _renderTestManagerPanel} from './package-testing';
import {functionSignatureEditor} from './function-signature-editor';
import {addToJSContextCommand, getMinifiedClassNameMap, _renderDevPanel} from './dev-panel';

import {_testDetectorsDialog, _testDetectorsStandard} from './utils/test-detectors';
import {testManagerViewNew} from './package-testing-new';

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

//name: renderTestManagerPanel
//tags: dev-tools
//input: object ent
//output: widget panel
export function renderTestManagerPanel(ent: EntityType): Promise<DG.Widget> {
  return _renderTestManagerPanel(ent);
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
    if (!acc.getPane('Test manager'))
      acc.addPane('Test manager', () => ui.wait(async () => (await renderTestManagerPanel(ent)).root));
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
  await testManagerView();
}

//name: TestManager1
//top-menu: Tools | Dev | Test Manager New
//tags: app
export async function testManager1(): Promise<void> {
  c = grok.functions.getCurrentCall();
  await testManagerViewNew();
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

//name: DummyFuncToInit
export function DummyFuncToInit() {
  return 0;
}
