import * as ui from 'datagrok-api/ui';
import * as Vue from 'vue';

const texts = new WeakMap<HTMLElement, string>();

// bound once as a getter: binding again cannot clear a tooltip
export const tooltip = {
  mounted: (el: HTMLElement, binding: Vue.DirectiveBinding<string>) => {
    texts.set(el, binding.value);
    ui.tooltip.bind(el, () => texts.get(el) ?? null);
  },
  updated: (el: HTMLElement, binding: Vue.DirectiveBinding<string>) => {
    texts.set(el, binding.value);
  },
};
