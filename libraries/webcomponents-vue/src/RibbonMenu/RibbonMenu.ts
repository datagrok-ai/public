import * as Vue from 'vue';
import {useViewService, RibbonMenuItem} from '../ViewService/ViewService';

export const RibbonMenu = Vue.defineComponent({
  name: 'RibbonMenu',
  props: {
    groupName: {
      type: String,
      required: true,
    },
    items: {
      type: Array as Vue.PropType<RibbonMenuItem[]>,
      required: true,
    },
    priority: {
      type: Number,
      required: false,
    },
  },
  setup(props) {
    const service = useViewService();
    service.registerMenuGroup(() => props.groupName, () => props.items, () => props.priority);
    return () => null;
  },
});
