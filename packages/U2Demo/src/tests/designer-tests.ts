import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {registerAll, renderSpec} from '@datagrok-libraries/u2';
import {SpecNodeRef, designerView, registerSpecNodeHandler} from '@datagrok-libraries/u2/src/dg/index.js';
import {DESIGNER_SPEC, designerContext} from '../designer';

/* The designer's platform seams — handler registration, selection, and the view's chrome. Everything
   here is client-side, so it runs under `grok test` and under the `?tests=` runner in local mode.
   Behaviour that needs a pointer (hit-test, adorners, mode toggle) is the e2e lane's job:
   libraries/u2/e2e. */
category('designer', () => {
  const designer = (): {ref: SpecNodeRef, dispose: () => void} => {
    registerAll();
    registerSpecNodeHandler();
    const instance = renderSpec(DESIGNER_SPEC, designerContext());
    const node = [...instance.nodes().keys()].find((n) => n.name === 'saveButton')!;
    return {ref: new SpecNodeRef(instance, node), dispose: () => instance.dispose()};
  };

  test('the node handler is registered and applicable', async () => {
    const {ref, dispose} = designer();
    try {
      const handler = DG.ObjectHandler.forEntity(ref);
      expect(handler != null, true);
      expect(handler!.type, 'u2-spec-node');
      expect(handler!.isApplicable(ref), true);
      expect(handler!.getCaption(ref), 'saveButton (u2-button)');
    } finally {
      dispose();
    }
  });

  // the context panel renders asynchronously and lags one assignment behind, so what is asserted is
  // the handler the platform resolves for the selected object — the same element the panel shows
  test('selection renders the caption', async () => {
    const {ref, dispose} = designer();
    try {
      grok.shell.o = ref;
      const properties = DG.ObjectHandler.forEntity(ref)!.renderProperties(ref);
      const caption = properties.querySelector('.u2-designer-caption');
      expect(caption != null, true);
      expect(caption!.textContent, 'saveButton (u2-button)');
      expect(properties.textContent!.includes('u2-button'), true);
    } finally {
      dispose();
    }
  });

  test('the view carries the structure tree and the status path', async () => {
    registerAll();
    const view = designerView(DESIGNER_SPEC, {name: 'Designer test', ctx: designerContext()});
    try {
      const labels = [...view.toolbox.querySelectorAll('.u2-tree-label')].map((l) => l.textContent);
      expect(labels.includes('layout'), true);
      const status = view.statusBarPanels.map((p) => p.textContent).join(' ');
      expect(/layout · \d+ nodes/.test(status), true);
    } finally {
      view.close();
    }
  });
});
