import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {useDebounceFn} from '@vueuse/core';

const LOADER_DEBOUNCE_TIME = 150;

export const ifOverlapping = {
  loaderMapping: new Map<HTMLElement, HTMLElement>(),
  updateFnMapping: new Map<HTMLElement, Function>(),

  mounted: (el: HTMLElement, binding: Vue.DirectiveBinding<boolean>) => {
    const customText = binding.arg;

    const loader = ui.divV([
      ui.label(customText ?? 'Updating...'),
      ui.loader(),
    ], 'd4-update-shadow');
    loader.style.zIndex = '1';
    ifOverlapping.loaderMapping.set(el, loader);

    const updateFn = (isOverlapping: boolean) => {
      const existingLoader = ifOverlapping.loaderMapping.get(el);
      if (isOverlapping && existingLoader) {
        el.append(existingLoader);
        el.classList.add('ui-box');
      }
      if (!isOverlapping && existingLoader) {
        existingLoader.remove();
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
  },
  updated: (el: HTMLElement, binding: Vue.DirectiveBinding<boolean>) => {
    const isOverlapping = binding.value;
    const debouncedFn = ifOverlapping.updateFnMapping.get(el)!;
    debouncedFn(isOverlapping);
  },
  beforeUnmount: (el: HTMLElement) => {
    const existingLoader = ifOverlapping.loaderMapping.get(el);
    if (existingLoader)
      ifOverlapping.loaderMapping.delete(el);

    const debouncedFn = ifOverlapping.updateFnMapping.get(el);
    if (debouncedFn)
      ifOverlapping.updateFnMapping.delete(el);
  },
};
