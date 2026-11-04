import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {createDgViewService} from '@datagrok-libraries/webcomponents-vue/src/ViewService/ViewService';
import {menuItem, panelItem} from './ribbon-test-utils';

category('WebComponentsVue: DG view service', () => {
  test('panels and menu groups reach the view and leave on dispose', async () => {
    const view = DG.View.create();
    const service = createDgViewService(view);
    const scope = Vue.effectScope();
    try {
      scope.run(() => {
        service.registerPanel(() => [panelItem({icon: 'save', tooltip: 'Save'})]);
        service.registerMenuGroup(() => 'Tools', () => [menuItem({text: 'Recalculate'})]);
      });
      await Vue.nextTick();
      const panels = view.getRibbonPanels();
      expect(panels.length, 1);
      expect(panels[0][0].querySelector('[aria-label="Save"]') != null, true);
      expect(!!view.ribbonMenu.find('Tools'), true, 'the menu group reached the view');

      service.dispose();
      expect(view.getRibbonPanels().length, 0);
      expect(!view.ribbonMenu.find('Tools'), true, 'the menu group stayed after dispose');
    } finally {
      scope.stop();
      view.close();
    }
  });
});
