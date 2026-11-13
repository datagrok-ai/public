import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {awaitCheck, category, expect, test} from '@datagrok-libraries/test/src/test';
import {Compute2ModelHandler} from '../model-handler';
import {BIOREACTOR_WORKFLOW, uniqueName, withRunStorage} from './run-storage';

category('Model Hub panel', () => {
  test('Shows the details, my last runs and the name to copy', () => withRunStorage(async (f) => {
    const workflow = DG.Func.byName(BIOREACTOR_WORKFLOW);
    const name = uniqueName('Model panel');
    await f.save(workflow.prepare({}), {metadata: {name}});
    const host = ui.div([new Compute2ModelHandler().renderProperties(workflow)]);
    document.body.append(host);
    try {
      const headers = () => [...host.querySelectorAll('.d4-accordion-pane-header')].map((h) => h.textContent?.trim());
      expect(JSON.stringify(headers()), JSON.stringify(['Details', 'My runs', 'Dev']));
      await awaitCheck(() => host.innerText.includes('Workflow') && host.innerText.includes(name),
        'Details or the saved run is not shown', 10000);
      expect(host.innerText.includes(workflow.nqName), false);
      (host.querySelectorAll('.d4-accordion-pane-header')[2] as HTMLElement).click();
      await awaitCheck(() => host.innerText.includes(workflow.nqName), 'The name is not shown in Dev', 5000);
    } finally {
      host.remove();
    }
  }));
});
