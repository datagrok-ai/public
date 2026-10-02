import * as Vue from 'vue';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {createRibbonService}
  from '@datagrok-libraries/webcomponents-vue/src/ViewService/ribbon-service';
import {countOp, fakeHost, identityDebounce, lastOp, menuItem, panelItem, runInScope} from './ribbon-test-utils';

function setup(initialPanels: HTMLElement[][] = []) {
  const {host, calls, setClosing} = fakeHost(initialPanels);
  const service = createRibbonService(host, {debounce: identityDebounce});
  return {service, calls, setClosing};
}

category('WebComponentsVue: Ribbon service', () => {
  test('writes registered panel once with external panels first', async () => {
    const external = [[document.createElement('div')]];
    const {service, calls} = setup(external);
    const items = Vue.ref([panelItem({icon: 'a'})]);
    const scope = runInScope(() => service.registerPanel(() => items.value));
    await Vue.nextTick();
    expectDeepEqual(countOp(calls, 'setPanels'), 1);
    const panels = lastOp(calls, 'setPanels').args[0] as HTMLElement[][];
    expectDeepEqual(panels.length, 2);
    expect(panels[0] === external[0]);
    expectDeepEqual(panels[1].length, 1);
    scope.stop();
  });

  test('batches same tick unregister and register into one write', async () => {
    const {service, calls} = setup();
    const scopeA = runInScope(() => service.registerPanel(() => [panelItem({icon: 'a'})]));
    await Vue.nextTick();
    expectDeepEqual(countOp(calls, 'setPanels'), 1);
    scopeA.stop();
    const scopeB = runInScope(() => service.registerPanel(() => [panelItem({icon: 'b'})]));
    await Vue.nextTick();
    expectDeepEqual(countOp(calls, 'setPanels'), 2);
    const panels = lastOp(calls, 'setPanels').args[0] as HTMLElement[][];
    expectDeepEqual(panels.length, 1);
    expect(panels[0][0].classList.contains('fa-b'));
    scopeB.stop();
  });

  test('skips host writes for property only changes', async () => {
    const {service, calls} = setup();
    const items = Vue.ref([panelItem({icon: 'x'})]);
    const scope = runInScope(() => service.registerPanel(() => items.value));
    await Vue.nextTick();
    const el = (lastOp(calls, 'setPanels').args[0] as HTMLElement[][])[0][0];
    items.value = [panelItem({icon: 'x', active: true})];
    await Vue.nextTick();
    expectDeepEqual(countOp(calls, 'setPanels'), 1);
    expectDeepEqual(el.style.backgroundColor, 'var(--grey-1)');
    items.value = [panelItem({icon: 'x'})];
    await Vue.nextTick();
    expectDeepEqual(countOp(calls, 'setPanels'), 1);
    expectDeepEqual(el.style.backgroundColor, '');
    scope.stop();
  });

  test('writes once per membership change and keeps element identity', async () => {
    const {service, calls} = setup();
    const items = Vue.ref([panelItem({icon: 'x'})]);
    const scope = runInScope(() => service.registerPanel(() => items.value));
    await Vue.nextTick();
    const el = (lastOp(calls, 'setPanels').args[0] as HTMLElement[][])[0][0];
    items.value = [panelItem({icon: 'x'}), panelItem({icon: 'y'})];
    await Vue.nextTick();
    expectDeepEqual(countOp(calls, 'setPanels'), 2);
    const panels = lastOp(calls, 'setPanels').args[0] as HTMLElement[][];
    expect(panels[0][0] === el);
    expectDeepEqual(panels[0].length, 2);
    scope.stop();
  });

  test('orders panels by priority with unset last', async () => {
    const {service, calls} = setup();
    const scopeA = runInScope(() => service.registerPanel(() => [panelItem({icon: 'a'})]));
    await Vue.nextTick();
    const scopeB = runInScope(() => service.registerPanel(() => [panelItem({icon: 'b'})], () => 0));
    await Vue.nextTick();
    const panels = lastOp(calls, 'setPanels').args[0] as HTMLElement[][];
    expect(panels[0][0].classList.contains('fa-b'));
    expect(panels[1][0].classList.contains('fa-a'));
    scopeA.stop();
    scopeB.stop();
  });

  test('merges same name menu groups across registrations', async () => {
    const {service, calls} = setup();
    const scope = runInScope(() => {
      service.registerMenuGroup(() => 'G', () => [menuItem({text: 'first'})]);
      service.registerMenuGroup(() => 'G', () => [menuItem({text: 'second'})]);
    });
    await Vue.nextTick();
    expectDeepEqual(countOp(calls, 'rebuildGroup'), 1);
    const [name, rank, els] = lastOp(calls, 'rebuildGroup').args;
    expectDeepEqual(name, 'G');
    expectDeepEqual(rank, 0);
    expectDeepEqual((els as HTMLElement[]).map((el) => el.textContent), ['first', 'second']);
    scope.stop();
  });

  test('removes group when last contributor leaves', async () => {
    const {service, calls} = setup();
    const scopeA = runInScope(() => service.registerMenuGroup(() => 'G', () => [menuItem({text: 'a'})]));
    const scopeB = runInScope(() => service.registerMenuGroup(() => 'G', () => [menuItem({text: 'b'})]));
    await Vue.nextTick();
    scopeB.stop();
    await Vue.nextTick();
    const els = lastOp(calls, 'rebuildGroup').args[2] as HTMLElement[];
    expectDeepEqual(els.map((el) => el.textContent), ['a']);
    expectDeepEqual(countOp(calls, 'removeGroup'), 0);
    scopeA.stop();
    await Vue.nextTick();
    expectDeepEqual(countOp(calls, 'removeGroup'), 1);
    expectDeepEqual(lastOp(calls, 'removeGroup').args[0], 'G');
  });

  test('moves items when group name changes', async () => {
    const {service, calls} = setup();
    const name = Vue.ref('G');
    const scope = runInScope(() => service.registerMenuGroup(() => name.value, () => [menuItem({text: 'a'})]));
    await Vue.nextTick();
    name.value = 'H';
    await Vue.nextTick();
    expectDeepEqual(lastOp(calls, 'removeGroup').args[0], 'G');
    expectDeepEqual(lastOp(calls, 'rebuildGroup').args[0], 'H');
    scope.stop();
  });

  test('serves check and disabled through pull callbacks without rebuild', async () => {
    const {service, calls} = setup();
    const checked = Vue.ref(false);
    const disabled = Vue.ref(false);
    const scope = runInScope(() => service.registerMenuGroup(() => 'G', () => [
      menuItem({text: 'a', check: checked.value, disabled: disabled.value, disabledReason: 'busy'}),
    ]));
    await Vue.nextTick();
    expectDeepEqual(countOp(calls, 'rebuildGroup'), 1);
    const [, , els, , opts] = lastOp(calls, 'rebuildGroup').args;
    const el = (els as HTMLElement[])[0];
    expect(opts.isChecked(el), false);
    expect(opts.isValid(el), null);
    checked.value = true;
    disabled.value = true;
    await Vue.nextTick();
    expectDeepEqual(countOp(calls, 'rebuildGroup'), 1);
    expect(opts.isChecked(el));
    expect(opts.isValid(el), 'busy');
    scope.stop();
  });

  test('dispatches menu clicks with disabled gating', async () => {
    const {service, calls} = setup();
    let ranFirst = 0;
    let ranSecond = 0;
    const scope = runInScope(() => service.registerMenuGroup(() => 'G', () => [
      menuItem({text: 'a', onClick: () => ranFirst++}),
      menuItem({text: 'b', onClick: () => ranSecond++,
        disabled: true, disabledReason: 'busy', disabledReasonMode: 'popup', disabledStyle: 'none'}),
    ]));
    await Vue.nextTick();
    const [, , els, dispatch] = lastOp(calls, 'rebuildGroup').args;
    dispatch(els[0]);
    expectDeepEqual(ranFirst, 1);
    dispatch(els[1]);
    expectDeepEqual(ranSecond, 0);
    expectDeepEqual(lastOp(calls, 'warn').args[0], 'busy');
    scope.stop();
  });

  test('suppresses writes while closing', async () => {
    const {service, calls, setClosing} = setup();
    const items = Vue.ref([panelItem({icon: 'a'})]);
    const scope = runInScope(() => service.registerPanel(() => items.value));
    await Vue.nextTick();
    expectDeepEqual(countOp(calls, 'setPanels'), 1);
    setClosing(true);
    items.value = [panelItem({icon: 'a'}), panelItem({icon: 'b'})];
    await Vue.nextTick();
    expectDeepEqual(countOp(calls, 'setPanels'), 1);
    scope.stop();
  });

  test('dispose restores external panels and removes groups', async () => {
    const external = [[document.createElement('div')]];
    const {service, calls} = setup(external);
    const scope = runInScope(() => {
      service.registerPanel(() => [panelItem({icon: 'a'})]);
      service.registerMenuGroup(() => 'G', () => [menuItem({text: 'a'})]);
    });
    await Vue.nextTick();
    service.dispose();
    const panels = lastOp(calls, 'setPanels').args[0] as HTMLElement[][];
    expectDeepEqual(panels.length, 1);
    expect(panels[0] === external[0]);
    expectDeepEqual(lastOp(calls, 'removeGroup').args[0], 'G');
    scope.stop();
  });

  test('stops reacting after dispose', async () => {
    const {service, calls} = setup();
    const items = Vue.ref([panelItem({icon: 'a'})]);
    const scope = runInScope(() => service.registerPanel(() => items.value));
    await Vue.nextTick();
    service.dispose();
    const before = calls.length;
    items.value = [panelItem({icon: 'b'})];
    await Vue.nextTick();
    expectDeepEqual(calls.length, before);
    scope.stop();
  });
});
