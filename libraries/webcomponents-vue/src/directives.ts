import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

export const ifOverlapping = {
  loaderMapping: new Map<HTMLElement, HTMLElement>(),

  mounted: (el: HTMLElement, binding: Vue.DirectiveBinding<boolean>) => {
    const customText = binding.arg;

    const loader = ui.divV([
      ui.label(customText ?? 'Updating...'),
      ui.loader(),
    ], 'd4-update-shadow');
    loader.style.zIndex = '1';
    ifOverlapping.loaderMapping.set(el, loader);
  },
  updated: (el: HTMLElement, binding: Vue.DirectiveBinding<boolean>) => {
    const isOverlapping = binding.value;

    const existingLoader = ifOverlapping.loaderMapping.get(el);
    if (isOverlapping && existingLoader) {
      el.append(existingLoader);
      el.classList.add('ui-box')
    }
    if (!isOverlapping && existingLoader) {
      existingLoader.remove();
      el.classList.remove('ui-box')
    }
  },
  beforeUnmount: (el: HTMLElement) => {
    const existingLoader = ifOverlapping.loaderMapping.get(el);
    if (existingLoader) {
      ifOverlapping.loaderMapping.delete(el);
    }
  },
};