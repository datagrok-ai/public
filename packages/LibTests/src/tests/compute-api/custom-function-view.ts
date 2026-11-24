import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {filter, take} from 'rxjs/operators';
import {awaitCheck, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {historyInput, openRunDialog, openWorkflowRun} from '@datagrok-libraries/compute-api';
import {historyUtils} from '@datagrok-libraries/compute-utils';
import {deleteRuns} from '../../test-utils';

// LibTests loads before Compute2, so the view class can only be declared once Compute2 is loaded
const makeApiView = () => class ApiView extends window.compute2.CustomFunctionView {
  loads = 0;

  buildIO() {
    return ui.div();
  }

  async onAfterLoadRun() {
    this.loads++;
  }
};

category('ComputeApi: custom function view', () => {
  before(async () => {
    await DG.Func.find({package: 'Compute2', name: 'init'})[0].prepare().call();
  });

  test('A view on Compute2\'s class saves and loads runs of its function', async () => {
    const ApiView = makeApiView();
    const view = new ApiView('LibTests:TestAdd2');
    await view.isReady.pipe(filter((ready) => ready), take(1)).toPromise();
    let saved: DG.FuncCall | undefined;
    try {
      expect(view.funcCall!.func.nqName, 'LibTests:TestAdd2');
      view.funcCall!.inputs['a'] = 2;
      view.funcCall!.inputs['b'] = 3;
      await view.funcCall!.call();
      saved = await view.saveRun(view.funcCall!, {name: `compute-api view ${Date.now()}`});
      expect(saved.outputs['res'], 5);
      await view.loadRun(saved.id);
      expect(view.loads, 1);
      expect(view.funcCall!.id, saved.id);
    } finally {
      await deleteRuns([saved?.id]);
    }
  });

  test('historyInput is Compute2\'s input', async () => {
    const input = historyInput('Previous run', {model: 'LibTests:TestAdd2'});
    expect(input instanceof DG.InputBase, true);
    expect(input.value, null);
    expect(input.caption, 'Previous run');
  });

  test('openRunDialog resolves to null when cancelled', async () => {
    const picked = openRunDialog(DG.Func.byName('LibTests:TestAdd2'));
    let browser: HTMLElement | undefined;
    await awaitCheck(() => !!(browser = document.querySelector<HTMLElement>('.d4-dialog .c2-runs-browser') ??
      undefined), 'the run browser should open', 15000);
    const dialog = browser!.closest<HTMLElement>('.d4-dialog')!;
    [...dialog.querySelectorAll('button')].find((b) => b.textContent?.toLowerCase().includes('cancel'))!.click();
    expect(await picked, null);
  });

  test('openWorkflowRun opens a saved run in its editor', async () => {
    const call = DG.Func.byName('LibTests:TestAdd2').prepare({a: 2, b: 3});
    await call.call();
    let saved: DG.FuncCall | undefined;
    let opened: DG.ViewBase | undefined;
    try {
      saved = await historyUtils.saveRun(call);
      await openWorkflowRun(saved.id);
      const isOpened = (v: DG.ViewBase) => v.path.includes(saved!.id);
      await awaitCheck(() => [...grok.shell.views].some(isOpened), 'the run should open', 15000);
      opened = [...grok.shell.views].find(isOpened);
    } finally {
      opened?.close();
      await deleteRuns([saved?.id]);
    }
  });
});
