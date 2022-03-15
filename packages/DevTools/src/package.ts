import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { scriptEditor } from "./script-editor";
import { IconTool } from "./icon-tool";
import { EntityType } from './constants';
import './styles.css';
import * as tests from "./tests/test-examples";
import {testPackages} from "./package-testing";
import {viewersGallery} from "./viewers-gallery";
import { functionSignatureEditor } from './function-signature-editor';
import { addToJSContextCommand, getMinifiedClassNameMap, _renderDevPanel } from './dev-panel';

export const _package = new DG.Package();
let minifiedClassNameMap = {};

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
    let devPane = acc.getPane('Dev');
    if (!devPane) devPane = acc.addPane('Dev', () => ui.wait(async () => (await renderDevPanel(ent)).root));
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

//description: ViewerGallery
//tags: autostart
export function _viewerGallery(): void { 
  //grok.shell.topMenu.find('Add').separator().item('Add viewer...', ()=>viewersGallery())
  grok.events.onViewAdded.subscribe((view) => {
    if (view.type == 'TableView'){
      let panel = view.getRibbonPanels();
      panel[0][1].remove();
      
      let icon = ui.iconFA('', ()=>{viewersGallery()});
      icon.className = 'grok-icon svg-icon svg-add-viewer';

      let btn = ui.div([icon]);
      btn.className = 'd4-ribbon-item';
      panel[0][0].after(btn)
    }
  });
}

//description: IconTool
export function _IconTool(): void { 
  grok.shell.newView('Icon Tool', [new IconTool('Icon Tool')]);
}

//name: testPackages
export async function _testPackages(): Promise<void> { await testPackages(); }

//tags: unitTest
export function _throwsException(): void { tests.throwsException(); }

//tags: unitTest
export function _returnsFine(): void { tests.returnsFine(); }