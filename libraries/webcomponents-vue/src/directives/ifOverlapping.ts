import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {useDebounceFn} from '@vueuse/core';

const LOADER_DEBOUNCE_TIME = 150;

type OverlappingState = {
  debouncedFn: (isOverlapping: boolean) => void,
  loader: HTMLElement,
  disposed: boolean,
};

export const ifOverlapping = {
  stateMapping: new Map<HTMLElement, OverlappingState>(),

  mounted: (el: HTMLElement, binding: Vue.DirectiveBinding<boolean>) => {
    const customText = binding.arg;

    const loader = ui.divV([
      ui.label(customText ?? 'Updating...'),
      ui.loader(),
    ], 'd4-update-shadow');
    loader.style.zIndex = '3000';

    const state: OverlappingState = {
      loader,
      disposed: false,
      // the disposed check keeps a show scheduled within the debounce window
      // from appending the loader to an already unmounted element
      debouncedFn: useDebounceFn((isOverlapping: boolean) => {
        if (state.disposed)
          return;
        if (isOverlapping)
          el.append(loader);
        else
          loader.remove();
      }, LOADER_DEBOUNCE_TIME),
    };
    ifOverlapping.stateMapping.set(el, state);

    ifOverlapping.updated(el, binding);
  },
  updated: (el: HTMLElement, binding: Vue.DirectiveBinding<boolean>) => {
    ifOverlapping.stateMapping.get(el)!.debouncedFn(binding.value);
  },
  beforeUnmount: (el: HTMLElement) => {
    const state = ifOverlapping.stateMapping.get(el);
    if (state) {
      state.disposed = true;
      state.loader.remove();
      ifOverlapping.stateMapping.delete(el);
    }
  },
};
