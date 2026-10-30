import * as Vue from 'vue';
import * as DG from 'datagrok-api/dg';
import {onClickOutside} from '@vueuse/core';
import {modelLabel} from './run-storage';

/** A searchable dropdown of models. */
export const ModelSelector = Vue.defineComponent({
  name: 'ModelSelector',
  props: {
    models: {type: Array as Vue.PropType<DG.Func[]>, required: true},
    modelValue: {type: DG.Func, required: false},
    isLoading: {type: Boolean, default: false},
    // an option for no model, shown first
    allLabel: {type: String, required: false},
  },
  emits: {'update:modelValue': (_func: DG.Func | undefined) => true},
  setup(props, {emit}) {
    const isOpen = Vue.ref(false);
    const search = Vue.ref('');
    const root = Vue.ref<HTMLElement>();
    onClickOutside(root, () => isOpen.value = false);

    const filtered = Vue.computed(() => {
      const query = search.value.toLowerCase();
      return props.models.filter((func) => !query || modelLabel(func).toLowerCase().includes(query));
    });

    const caption = () => props.modelValue ? modelLabel(props.modelValue) :
      (props.allLabel ?? (props.isLoading ? 'Loading models…' : 'Select model...'));

    return () => (
      <div ref={root} class='c2-runs-model-selector'>
        <div class='c2-runs-model-current' onClick={() => {isOpen.value = !isOpen.value; search.value = '';}}>
          <span>{caption()}</span>
          <span>&#9662;</span>
        </div>
        { isOpen.value &&
          <div class='c2-runs-model-dropdown'>
            <input type='text' class='c2-runs-model-search' placeholder='Filter models...' value={search.value}
              onInput={(e: Event) => search.value = (e.target as HTMLInputElement).value}/>
            { props.allLabel && !search.value &&
              <div class='c2-runs-row' onClick={() => {emit('update:modelValue', undefined); isOpen.value = false;}}>
                {props.allLabel}
              </div> }
            { filtered.value.map((func) =>
              <div key={func.id} class='c2-runs-row'
                onClick={() => {emit('update:modelValue', Vue.markRaw(func)); isOpen.value = false;}}>
                {modelLabel(func)}
              </div>) }
            { filtered.value.length === 0 && <div class='c2-runs-muted c2-runs-pad'>No models found</div> }
          </div> }
      </div>
    );
  },
});
