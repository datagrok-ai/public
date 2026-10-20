import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {category, test, expect, awaitCheck} from '@datagrok-libraries/test/src/test';
import {injectInputBaseStatus} from '@datagrok-libraries/webcomponents-vue/src/InputForm/utils';
import {awaitWebComponents} from './utils';

const ADVICE = {validation: {warnings: [{
  description: 'Too high',
  actions: [{actionName: 'Set to 10', action: 'action-1', additionalParams: {val: 10}}],
}]}};

async function clickAdviceLink(root: HTMLElement) {
  const icon = () => root.querySelector('dg-validation-icon i') as HTMLElement | null;
  await awaitCheck(() => icon() != null, 'validation icon not rendered', 3000);
  icon()!.click();
  const link = () => Array.from(root.querySelectorAll('a.ui-link'))
    .find((el) => el.textContent?.trim() === 'Set to 10') as HTMLElement | undefined;
  await awaitCheck(() => link() != null, 'advice link not rendered', 3000);
  link()!.click();
}

category('Validation advice actions', () => {
  test('Advice link sends its additional params with the action request', async () => {
    await awaitWebComponents();
    const host = ui.div();
    document.body.append(host);
    try {
      const icon = document.createElement('dg-validation-icon') as any;
      host.append(icon);
      const requests: any[] = [];
      icon.addEventListener('action-request', (ev: any) => requests.push([ev.detail, ev.additionalParams]));
      icon.validationStatus = ADVICE;
      await clickAdviceLink(host);
      expect(JSON.stringify(requests), JSON.stringify([['action-1', {val: 10}]]));
    } finally {
      host.remove();
    }
  });

  test('Input form forwards advice params with the action request', async () => {
    await awaitWebComponents();
    const input = ui.input.forProperty(DG.Property.js('x', DG.TYPE.FLOAT));
    const emitted: any[] = [];
    injectInputBaseStatus((...args: any[]) => emitted.push(args), 'x', input);
    document.body.append(input.root);
    try {
      (input as any).setStatus(ADVICE);
      await clickAdviceLink(input.root);
      expect(JSON.stringify(emitted), JSON.stringify([['actionRequested', 'action-1', {val: 10}]]));
    } finally {
      input.root.remove();
    }
  });
});
