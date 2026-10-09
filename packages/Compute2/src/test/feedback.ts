import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {category, test} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {useFeedbackItems} from '../composables/use-feedback';

const withSettings = (settings?: Record<string, string>) => ({package: {settings}}) as unknown as DG.Func;

category('Feedback menu', () => {
  test('items come from the package feedback urls', async () => {
    const func = Vue.shallowRef<DG.Func | undefined>(
      withSettings({REPORT_BUG_URL: 'https://bugs', REQUEST_FEATURE_URL: 'https://features'}));
    const items = useFeedbackItems(func);
    expectDeepEqual(items.value.map((item) => item.text), ['Report a bug', 'Request a feature']);
    const opened: string[] = [];
    const open = window.open;
    window.open = ((url: string) => {
      opened.push(url);
      return null;
    }) as typeof window.open;
    try {
      items.value.forEach((item) => item.onClick?.());
    } finally {
      window.open = open;
    }
    expectDeepEqual(opened, ['https://bugs', 'https://features']);

    func.value = withSettings({REQUEST_FEATURE_URL: 'https://features'});
    expectDeepEqual(items.value.map((item) => item.text), ['Request a feature'], {prefix: 'One url'});
    func.value = withSettings();
    expectDeepEqual(items.value, [], {prefix: 'No settings'});
    func.value = undefined;
    expectDeepEqual(items.value, [], {prefix: 'No function'});
  });
});
