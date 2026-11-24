import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {filter, take} from 'rxjs/operators';
import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {historyInput} from '@datagrok-libraries/compute-api';
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
});
