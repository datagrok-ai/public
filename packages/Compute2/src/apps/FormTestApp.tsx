/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {InputForm} from '@datagrok-libraries/webcomponents-vue';
import * as Vue from 'vue';

export const FormTestApp = Vue.defineComponent({
  name: 'FormTestApp',
  setup() {
    const fc = Vue.shallowRef<DG.FuncCall | undefined>(undefined);

    const logFuncCall = () => {
      console.log(Object.entries(fc.value?.inputs ?? {}));
    };

    const changeFuncCall = () => {
      const nfc = fc.value?.func.nqName === 'Compute2:SimpleInputs2' ?
        DG.Func.byName('Compute2:ObjectCooling2').prepare():
        DG.Func.byName('Compute2:SimpleInputs2').prepare();
      fc.value = Vue.markRaw(nfc);
      Vue.triggerRef(fc);
    };

    const removeFuncCall = () => {
      fc.value = undefined;
    };

    return () => (
      <div style={{width: '100%', height: '100%'}}>
        <button onClick={logFuncCall}>log funccall inputs</button>
        <button onClick={changeFuncCall}>change funcall</button>
        <button onClick={removeFuncCall}>remove funcall</button>
        { fc.value && <InputForm funcCall={fc.value}></InputForm> }
      </div>
    );
  },
});
