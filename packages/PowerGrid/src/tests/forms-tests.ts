import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {waitUntil, category, expect, test} from '@datagrok-libraries/test/src/test';
import type {FormsViewer} from '@datagrok-libraries/utils/src/viewers/forms-viewer';
import {layoutRoundTrip, projectRoundTrip} from './utils';

const FIELDS = ['age', 'subj', 'height', 'sex'];
const PINNED_ROW = 7;

/** The viewer's own object: the platform hands out a plain viewer until the package has created it. */
async function formsViewer(view: DG.TableView): Promise<FormsViewer> {
  const find = () =>
    Array.from(view.viewers).find((v) => typeof (v as any).pinRow === 'function') as unknown as FormsViewer;
  await waitUntil(() => find() != null, 'no Forms viewer in the view', 5000);
  return find();
}

async function formsView(): Promise<DG.TableView> {
  const view = grok.shell.addTableView(grok.data.demo.demog(100));
  view.addViewer('Forms');
  const forms = await formsViewer(view);
  forms.setOptions({fieldsColumnNames: FIELDS, sortBy: 'age'});
  forms.pinRow(PINNED_ROW, 'subj');
  forms.render();
  return view;
}

async function expectForms(view: DG.TableView): Promise<void> {
  const forms = await formsViewer(view);
  const pinned = String(view.dataFrame.get('subj', PINNED_ROW));
  const values = () => forms.getWidgetStatus().values!;
  await waitUntil(() => !forms.isRenderPending && values()['pinned records'] === 1, 'no pinned card', 5000);
  expect(values()['fields'], FIELDS.join(', '));
  expect(values()['sort column'], 'age');
  expect(values()['subj of pinned card 1'], pinned);
  expect(forms.pinnedRowValues.join(), pinned);
  expect(forms.pinnedRowIndexes.join(), String(PINNED_ROW));
}

category('Forms viewer', () => {
  test('Layout', async () => {
    await layoutRoundTrip(await formsView(), grok.data.demo.demog(100), expectForms);
  });

  test('Project', async () => {
    await projectRoundTrip(await formsView(), expectForms);
  });

  test('A pinned row stays pinned under a filter that hides it and leaves the pinned pane once deleted', async () => {
    const view = await formsView();
    try {
      const forms = await formsViewer(view);
      const df = view.dataFrame;
      const pinned = () => forms.getWidgetStatus().values!['pinned records'];
      await waitUntil(() => !forms.isRenderPending && pinned() === 1, 'no pinned card', 5000);
      df.filter.init((i) => i !== PINNED_ROW);
      await waitUntil(() => !forms.isRenderPending, 'the cards did not settle after the filter', 5000);
      expect(pinned(), 1, 'a pinned row the filter hides left the pinned pane');
      df.filter.setAll(true);
      df.rows.removeAt(PINNED_ROW);
      await waitUntil(() => !forms.isRenderPending && pinned() === 0, 'a deleted pinned row stayed in the pinned pane', 5000);
      expect(forms.pinnedRowValues.length, 0);
    }
    finally {
      grok.shell.closeAll();
    }
  });
});
