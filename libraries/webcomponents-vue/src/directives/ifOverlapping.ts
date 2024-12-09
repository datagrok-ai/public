import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {useDebounceFn} from '@vueuse/core';

const LOADER_DEBOUNCE_TIME = 150;

export const ifOverlapping = {
  updateFnMapping: new Map<HTMLElement, Function>(),

  mounted: (el: HTMLElement, binding: Vue.DirectiveBinding<boolean>) => {
    const customText = binding.arg;

    const loader = ui.divV([
      ui.label(customText ?? 'Updating...'),
      ui.loader(),
    ], 'd4-update-shadow');
    loader.style.zIndex = '1';

    const updateFn = (isOverlapping: boolean) => {
      if (isOverlapping && loader) {
        el.append(loader);
        el.classList.add('ui-box');
      }
      if (!isOverlapping && loader) {
        loader.remove();
        el.classList.remove('ui-box');
      }
    };

    ifOverlapping.updateFnMapping.set(
      el,
      useDebounceFn(
        (isOverlapping: boolean) => updateFn(isOverlapping),
        LOADER_DEBOUNCE_TIME,
      ),
    );

    ifOverlapping.updated(el, binding);
  },
  updated: (el: HTMLElement, binding: Vue.DirectiveBinding<boolean>) => {
    const isOverlapping = binding.value;
    const debouncedFn = ifOverlapping.updateFnMapping.get(el)!;
    debouncedFn(isOverlapping);
  },
  beforeUnmount: (el: HTMLElement) => {
    const debouncedFn = ifOverlapping.updateFnMapping.get(el);
    if (debouncedFn)
      ifOverlapping.updateFnMapping.delete(el);
  },
};
