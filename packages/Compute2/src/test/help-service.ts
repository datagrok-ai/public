import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {awaitCheck, category, expect, test} from '@datagrok-libraries/test/src/test';
import {useHelp} from '../composables/use-help';

category('Help service', () => {
  test('Help keeps working after a load that fails', async () => {
    const scope = Vue.effectScope();
    const help = scope.run(() => useHelp())!;
    try {
      // a func whose id can't be read makes the help load fail
      const failing = {options: {help: 'x.md'}, get id(): string {
        throw new Error('no id');
      }} as unknown as DG.Func;
      help.changeHelpFunc(failing);
      await awaitCheck(() => !help.helpLoading.value, 'the failed load should finish', 3000);
      expect(help.helpContent.value == null, true);

      // the load is asynchronous, so a stream still running shows it as loading right away
      help.changeHelpFunc(failing);
      expect(help.helpLoading.value, true);
      await awaitCheck(() => !help.helpLoading.value, 'a later load should finish', 3000);
    } finally {
      scope.stop();
    }
  });
});
