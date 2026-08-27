/** Home-screen template cards + the bundled Interactive Viewers demo.
 *  The card click reads `files/` (AppData); the graph itself is ALSO verified
 *  through the deployed script entity (DB-backed), so a stand whose AppData
 *  file serving is down still proves the flow deploys and fully deserializes. */
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {FuncFlowView} from '../funcflow-view';
import {FlowEditor} from '../rete/flow-editor';
import {registerBuiltinNodes, registerAllFunctions} from '../rete/node-factory';
import {parseFlowBody, FLOW_LANGUAGE} from '../serialization/flow-script-format';
import {deserializeFlow} from '../serialization/flow-serializer';
import {_package} from '../package';
import {makeEditor, destroyEditor, until} from './test-utils';

const TEMPLATE_SLUGS = ['workflow-demo', 'sequence-demo', 'interactive-viewers'];
const IV_NODES = 16;
const IV_CONNECTIONS = 13;

category('Flow: start templates', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('every bundled template has a card on the home screen', async () => {
    const view = new FuncFlowView();
    try {
      await until(() => (view as any).flow != null, 10000);
      for (const slug of TEMPLATE_SLUGS)
        expect(!!view.root.querySelector(`[data-testid="ff-start-template-${slug}"]`), true, `${slug} card shown`);
    } finally {
      view.detach();
      view.root.remove();
    }
  });

  test('the deployed Interactive Viewers script parses and fully deserializes', async () => {
    const flows = await grok.dapi.scripts.filter(`language = "${FLOW_LANGUAGE}"`).list() as DG.Script[];
    const script = flows.find((s) => (s.friendlyName === 'Interactive Viewers' || s.name === 'InteractiveViewers'));
    expect(script != null, true,
      `the bundled scripts/Interactive Viewers.flow deployed as a script entity; saw: ${
        flows.map((s) => s.name).join(', ')}`);
    const {header, doc} = parseFlowBody(script!.script);
    expect(header.includes('//name: Interactive Viewers'), true, 'canonical name line');
    expect(header.includes('//input: string molecule'), true, 'the sketcher input surfaces in the header');
    expect(doc.metadata?.settings?.autorun, true, 'saved live — reopens live');
    const e = makeEditor();
    try {
      await deserializeFlow(doc, e.flow);
      expect(e.flow.getNodes().length, IV_NODES,
        'every node type resolved (an unregistered func is silently skipped)');
      expect(e.flow.getConnections().length, IV_CONNECTIONS, 'all wires survived');
    } finally {
      destroyEditor(e);
    }
  }, {timeout: 30000});

  test('the card click loads the whole flow (needs AppData file serving)', async () => {
    try {
      await _package.files.readAsText('Workflow Demo.flow');
    } catch {
      console.warn('Flow: start templates: AppData file serving is down on this stand ' +
        '(the pre-existing Workflow Demo.flow is unreadable too) — card-click load skipped');
      expect(true, true, 'skipped: AppData files unavailable');
      return;
    }
    const view = new FuncFlowView();
    try {
      await until(() => (view as any).flow != null, 10000);
      view.root.querySelector<HTMLElement>('[data-testid="ff-start-template-interactive-viewers"]')!.click();
      const flow = (view as any).flow as FlowEditor;
      expect(await until(() => flow.getNodes().length === IV_NODES, 20000), true,
        `template loaded — got ${flow.getNodes().length}/${IV_NODES} nodes`);
      expect(flow.getConnections().length, IV_CONNECTIONS, 'all wires survived');
      const overlay = view.root.querySelector<HTMLElement>('[data-testid="ff-start-overlay"]');
      expect(overlay?.style.display, 'none', 'the start panel is gone after loading');
      expect(((view as any).autorunScheduler as {enabled: boolean}).enabled, true,
        'the template reopens with autorun on');
    } finally {
      view.detach();
      view.root.remove();
    }
  }, {timeout: 60000});
});
