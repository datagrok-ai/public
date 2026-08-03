import * as Vue from 'vue';
import {RunComparison} from '../components/RunComparison/RunComparison';

export const RunComparisonApp = Vue.defineComponent({
  name: 'RunComparisonApp',
  props: {
    roleOnlyFilter: {
      type: Boolean,
      default: false,
    },
  },
  setup(props) {
    return () => <RunComparison roleOnlyFilter={props.roleOnlyFilter}/>;
  },
});
