import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {filter, take} from 'rxjs/operators';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {CustomFunctionView} from '@datagrok-libraries/compute-api';
import {deleteRuns} from '../../test-utils';

class ApiView extends CustomFunctionView {
  loads = 0;

  buildIO() {
    return ui.div();
  }

  async onAfterLoadRun() {
    this.loads++;
  }
}

category('ComputeApi: custom function view', () => {
  test('A view on the exported class saves and loads runs of its function', async () => {
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
});
