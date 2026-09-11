import * as Vue from 'vue';
import {useViewService, RibbonPanelItem} from '../ViewService/ViewService';

export const RibbonPanel = Vue.defineComponent({
  name: 'RibbonPanel',
  props: {
    items: {
      type: Array as Vue.PropType<RibbonPanelItem[]>,
      required: true,
    },
    priority: {
      type: Number,
      required: false,
    },
  },
  setup(props) {
    const service = useViewService();
    service.registerPanel(() => props.items, () => props.priority);
    return () => null;
  },
});
