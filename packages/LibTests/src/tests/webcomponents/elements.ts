import {category, test, expect, waitUntil, delay} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import type {DGIconFA, DGToggleInput} from '@datagrok-libraries/webcomponents/src/elements';
import {awaitWebComponents} from '../../test-utils';

async function mount<T extends HTMLElement>(tag: string) {
  await awaitWebComponents();
  const el = document.createElement(tag) as T;
  document.body.append(el);
  return {el, dispose: () => el.remove()};
}

category('WebComponents: elements', () => {
  test('the icon applies its name, style, cursor and animation', async () => {
    const {el, dispose} = await mount<DGIconFA>('dg-icon-fa');
    try {
      el.name = 'cog';
      el.faStyle = 'fas';
      el.animation = 'spin';
      el.cursor = 'move';
      const icon = el.querySelector('i')!;
      expect(icon.classList.contains('fa-cog'), true);
      expect(icon.classList.contains('fas') && !icon.classList.contains('fal'), true);
      expect(icon.classList.contains('fa-spin'), true);
      expect(icon.style.cursor, 'move');
      expect(el.childElementCount, 1);
    } finally {
      dispose();
    }
  });

  test('the toggle announces user changes', async () => {
    const {el, dispose} = await mount<DGToggleInput>('dg-toggle-input');
    try {
      const values: boolean[] = [];
      el.caption = 'Shared scale';
      el.addEventListener('value-changed', (ev: any) => values.push(ev.detail));
      el.querySelector<HTMLElement>('.ui-input-switch')!.click();
      await waitUntil(() => values.length > 0, 'a click did not announce the change', 3000);
      expectDeepEqual(values, [true]);
    } finally {
      dispose();
    }
  });

  test('a value set by the host is not announced back', async () => {
    const {el, dispose} = await mount<DGToggleInput>('dg-toggle-input');
    try {
      const values: boolean[] = [];
      el.caption = 'Shared scale';
      el.addEventListener('value-changed', (ev: any) => values.push(ev.detail));
      el.value = true;
      await delay(50);
      expectDeepEqual(values, [], {prefix: 'the host value came back as a user change'});
    } finally {
      dispose();
    }
  });
});
