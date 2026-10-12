import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {awaitCheck, category, expect, test} from '@datagrok-libraries/test/src/test';
import type {FormsViewer} from '@datagrok-libraries/utils/src/viewers/forms-viewer';
import {layoutRoundTrip, projectRoundTrip} from './utils';

const FIELDS = ['age', 'subj', 'height', 'sex'];
const PINNED_ROW = 7;

/** The viewer's own object: the platform hands out a plain viewer until the package has created it. */
async function formsViewer(view: DG.TableView): Promise<FormsViewer> {
  const find = () =>
    Array.from(view.viewers).find((v) => typeof (v as any).pinRow === 'function') as unknown as FormsViewer;
  await awaitCheck(() => find() != null, 'no Forms viewer in the view', 5000);
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
  await awaitCheck(() => !forms.isRenderPending && values()['pinned records'] === 1, 'no pinned card', 5000);
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
});
