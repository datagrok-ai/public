import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

export const ifOverlapping = {
  loaderMapping: new Map<HTMLElement, HTMLElement>(),

  updated: (el: HTMLElement, binding: Vue.DirectiveBinding<boolean>) => {
    const isOverlapping = binding.value;
    const customText = binding.arg;
    const existingLoader = ifOverlapping.loaderMapping.get(el);
    if (isOverlapping && !existingLoader) {
      const loader = ui.divV([
        ui.label(customText ?? 'Updating...'),
        ui.loader(),
      ], 'd4-update-shadow');
      el.append(loader);
      el.classList.add('ui-box')
      ifOverlapping.loaderMapping.set(el, loader);
    }
    if (!isOverlapping && existingLoader) {
      ifOverlapping.loaderMapping.delete(el);
      existingLoader.remove();
      el.classList.remove('ui-box')
    }
  },
};