import * as Vue from 'vue';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {createRibbonService}
  from '@datagrok-libraries/webcomponents-vue/src/ViewService/ribbon-service';
import {DgViewService, DG_VIEW_SERVICE_KEY, useViewService}
  from '@datagrok-libraries/webcomponents-vue/src/ViewService/ViewService';
import {RibbonPanel} from '@datagrok-libraries/webcomponents-vue/src/RibbonPanel/RibbonPanel';
import {RibbonMenu} from '@datagrok-libraries/webcomponents-vue/src/RibbonMenu/RibbonMenu';
import {countOp, fakeHost, identityDebounce, lastOp, menuItem, panelItem} from './ribbon-test-utils';

function setup() {
  const {host, calls} = fakeHost();
  const service = createRibbonService(host, {debounce: identityDebounce});
  return {service: {...service, view: undefined as any} as DgViewService, calls};
}

function mount(render: () => any, service: DgViewService) {
  const app = Vue.createApp({render});
  app.config.warnHandler = () => {};
  app.provide(DG_VIEW_SERVICE_KEY, service);
  app.mount(document.createElement('div'));
  return app;
}

category('WebComponentsVue: Ribbon components', () => {
  test('panel component registers and unregisters through the service', async () => {
    const {service, calls} = setup();
    const app = mount(() => Vue.h(RibbonPanel, {items: [panelItem({icon: 'a'})]}), service);
    await Vue.nextTick();
    expectDeepEqual(countOp(calls, 'setPanels'), 1);
    expectDeepEqual((lastOp(calls, 'setPanels').args[0] as HTMLElement[][]).length, 1);
    app.unmount();
    await Vue.nextTick();
    expectDeepEqual(countOp(calls, 'setPanels'), 2);
    expectDeepEqual((lastOp(calls, 'setPanels').args[0] as HTMLElement[][]).length, 0);
  });

  test('menu component registers its group', async () => {
    const {service, calls} = setup();
    const app = mount(() => Vue.h(RibbonMenu, {groupName: 'G', items: [menuItem({text: 'a'})]}), service);
    await Vue.nextTick();
    expectDeepEqual(countOp(calls, 'rebuildGroup'), 1);
    expectDeepEqual(lastOp(calls, 'rebuildGroup').args[0], 'G');
    app.unmount();
    await Vue.nextTick();
    expectDeepEqual(lastOp(calls, 'removeGroup').args[0], 'G');
  });

  test('items prop replacement flows through', async () => {
    const {service, calls} = setup();
    const items = Vue.ref([panelItem({icon: 'a'})]);
    const app = mount(() => Vue.h(RibbonPanel, {items: items.value}), service);
    await Vue.nextTick();
    expectDeepEqual(countOp(calls, 'setPanels'), 1);
    items.value = [panelItem({icon: 'a'}), panelItem({icon: 'b'})];
    await Vue.nextTick();
    expectDeepEqual(countOp(calls, 'setPanels'), 2);
    expectDeepEqual((lastOp(calls, 'setPanels').args[0] as HTMLElement[][])[0].length, 2);
    app.unmount();
  });

  test('useViewService throws without a provider', async () => {
    let error: unknown;
    const app = Vue.createApp({
      setup() {
        try {
          useViewService();
        } catch (e) {
          error = e;
        }
        return () => null;
      },
    });
    app.config.warnHandler = () => {};
    app.mount(document.createElement('div'));
    app.unmount();
    expect(error instanceof Error);
    expect((error as Error).message.includes('provideDgViewService'));
  });
});
