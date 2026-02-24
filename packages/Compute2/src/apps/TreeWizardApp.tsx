import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import '@he-tree/vue/style/default.css';
import '@he-tree/vue/style/material-design.css';
import {TreeWizard} from '../components/TreeWizard/TreeWizard';
import {PipelineInstanceConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';

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
    view: {
      type: DG.View,
      required: true,
    },
    resolve: {
      type: Function,
      required: false,
    },
  },
  setup(props) {
    const currentView = Vue.computed(() => Vue.markRaw(props.view));
    const resolve = Vue.computed(() => props.resolve ? Vue.markRaw(props.resolve) : undefined);
    const onReturn = (data: any) => {
      if (resolve.value)
        resolve.value(data);
      if (currentView.value)
        currentView.value.close();
    };
    return () => (
      <TreeWizard providerFunc={props.providerFunc} version={props.version} instanceConfig={props.instanceConfig} modelName={props.modelName} view={currentView.value} showReturn={!!resolve.value} onReturn={onReturn}/>
    );
  },
});
