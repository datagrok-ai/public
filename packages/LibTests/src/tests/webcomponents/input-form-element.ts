import * as DG from 'datagrok-api/dg';
import {category, test, expect, waitUntil, delay} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import type {InputForm} from '@datagrok-libraries/webcomponents/src/InputForm/InputForm';
import {awaitWebComponents} from '../../test-utils';
import {scriptWith, typeInto} from './fixtures';

async function mountForm() {
  await awaitWebComponents();
  const el = document.createElement('dg-input-form') as InputForm;
  const replaced: (DG.InputForm | undefined)[] = [];
  const changes: unknown[] = [];
  const validity: boolean[] = [];
  el.addEventListener('form-replaced', (ev: any) => replaced.push(ev.detail));
  el.addEventListener('input-changed', (ev: any) => changes.push(ev.detail));
  el.addEventListener('validation-changed', (ev: any) => validity.push(ev.detail));
  document.body.append(el);
  const dispose = () => {
    el.destroy();
    el.remove();
  };
  return {el, replaced, changes, validity, dispose};
}

category('WebComponents: InputForm element', () => {
  test('rapid funcCall swaps leave only the latest form', async () => {
    const {el, replaced, dispose} = await mountForm();
    try {
      const first = scriptWith('double a').prepare({a: 1});
      const second = scriptWith('double b').prepare({b: 2});
      el.funcCall = first;
      el.funcCall = second;
      await waitUntil(() => replaced.length > 0, 'form not built', 5000);
      await delay(300);
      expect(replaced.length, 1, 'only the latest form is announced');
      expect(el.funcCall === second, true);
      expect(el.childElementCount, 1);
      expect(!!replaced[0]!.getInput('b'), true);
    } finally {
      dispose();
    }
  });

  test('a cleared funcCall keeps the event streams alive', async () => {
    const {el, replaced, changes, dispose} = await mountForm();
    try {
      el.funcCall = scriptWith('double a').prepare({a: 1});
      await waitUntil(() => replaced.length === 1, 'first form not built', 5000);
      el.funcCall = undefined;
      await waitUntil(() => replaced.length === 2, 'clearing not announced', 5000);
      expect(replaced[1] == null, true);
      expect(el.childElementCount, 0);
      el.funcCall = scriptWith('double b').prepare({b: 2});
      await waitUntil(() => replaced.length === 3, 'second form not built', 5000);
      typeInto(replaced[2]!.getInput('b'), '5');
      await waitUntil(() => changes.length > 0, 'input changes no longer reach the element', 3000);
    } finally {
      dispose();
    }
  });

  test('a form opened with a missing required value reports itself invalid', async () => {
    const {el, validity, dispose} = await mountForm();
    try {
      el.funcCall = scriptWith('string s').prepare();
      await waitUntil(() => validity.length > 0, 'validity not reported', 5000);
      expectDeepEqual(validity, [false]);
      el.funcCall = scriptWith('string s').prepare({s: 'x'});
      await waitUntil(() => validity.length > 1, 'validity of the second form not reported', 5000);
      expectDeepEqual(validity, [false, true]);
    } finally {
      dispose();
    }
  });

  test('destroy cancels a form still being built', async () => {
    const {el, replaced, dispose} = await mountForm();
    try {
      el.funcCall = scriptWith('double a').prepare({a: 1});
      el.destroy();
      await delay(500);
      expect(el.childElementCount, 0, 'the pending form landed in the destroyed element');
      expect(replaced.length, 0, 'the destroyed element announced the pending form');
    } finally {
      dispose();
    }
  });
});
