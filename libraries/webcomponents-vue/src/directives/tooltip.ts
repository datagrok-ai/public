import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

export const tooltip = {
  mounted: (el: HTMLElement, binding: Vue.DirectiveBinding<string>) => {
    const tooltipText = binding.value;

    ui.tooltip.bind(el, tooltipText);
  },
};
