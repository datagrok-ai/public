import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import '@he-tree/vue/style/default.css';
import '@he-tree/vue/style/material-design.css';
import {TreeWizard} from '../components/TreeWizard/TreeWizard';
import {PipelineInstanceConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {useDgView} from '@datagrok-libraries/webcomponents-vue';

export const TreeWizardApp = Vue.defineComponent({
  name: 'TreeWizardApp',
  props: {
    providerFunc: {
      type: String,
      required: true,
    },
    version: {
      type: String,
    },
    modelName: {
      type: String,
      required: true,
    },
    instanceConfig: {
      type: Object as Vue.PropType<PipelineInstanceConfig>,
      required: false,
    },
    initialRunId: {
      type: String,
      required: false,
    },
    resolve: {
      type: Function,
      required: false,
    },
  },
  setup(props) {
    const currentView = useDgView();
    const resolve = Vue.computed(() => props.resolve ? Vue.markRaw(props.resolve) : undefined);
    const onReturn = (data: any) => {
      if (resolve.value)
        resolve.value(data);
      currentView.close();
    };
    return () => (
      <TreeWizard providerFunc={props.providerFunc} version={props.version} instanceConfig={props.instanceConfig} initialRunId={props.initialRunId} modelName={props.modelName} showReturn={!!resolve.value} onReturn={onReturn}/>
    );
  },
});
