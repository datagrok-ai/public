import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

export const ifOverlapping = {
  loaderMapping: new Map<HTMLElement, HTMLElement>(),

  updated: (el: HTMLElement, binding: Vue.DirectiveBinding<boolean>) => {
    const isOverlapping = binding.value;
    const existingLoader = ifOverlapping.loaderMapping.get(el);
    if (isOverlapping && !existingLoader) {
      const loader = ui.divV([
        ui.label('Updating...'),
        ui.loader(),
      ], 'd4-update-shadow');
      el.append(loader);
      ifOverlapping.loaderMapping.set(el, loader);
    }
    if (!isOverlapping && existingLoader) 
      existingLoader.remove();
    
  },
};