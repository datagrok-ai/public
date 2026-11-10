import * as Vue from 'vue';
import {tabWithin} from '../../composables/use-keyboard';

/**
 * A run's name, description and tags. Tab goes round the fields wherever they are shown, Enter in a one-line field
 * submits, and Name takes the focus when they open.
 */
export const RunMetadataFields = Vue.defineComponent({
  name: 'RunMetadataFields',
  props: {
    name: {type: String, required: true},
    description: {type: String, default: ''},
    tags: {type: String, default: ''},
    showDetails: {type: Boolean, default: true},
  },
  emits: {
    'update:name': (_value: string) => true,
    'update:description': (_value: string) => true,
    'update:tags': (_value: string) => true,
    'submit': () => true,
  },
  setup(props, {emit}) {
    const root = Vue.ref<HTMLElement>();
    const nameInput = Vue.ref<HTMLInputElement>();
    Vue.onMounted(() => nameInput.value?.focus());
    const onEnter = (e: KeyboardEvent) => e.key === 'Enter' && emit('submit');

    return () => (
      <div class='c2-runs-form' ref={root} onKeydown={(e: KeyboardEvent) => tabWithin(e, [root.value!])}>
        <label class='c2-runs-field'><span>Name</span>
          <input class='ui-input-editor' ref={nameInput} value={props.name} onKeydown={onEnter}
            onInput={(e: Event) => emit('update:name', (e.target as HTMLInputElement).value)}/></label>
        { props.showDetails && <label class='c2-runs-field'><span>Description</span>
          <textarea class='ui-input-editor' rows={5} value={props.description}
            onInput={(e: Event) => emit('update:description', (e.target as HTMLTextAreaElement).value)}/></label> }
        { props.showDetails && <label class='c2-runs-field'><span>Tags</span>
          <input class='ui-input-editor' placeholder='Comma-separated' value={props.tags} onKeydown={onEnter}
            onInput={(e: Event) => emit('update:tags', (e.target as HTMLInputElement).value)}/></label> }
      </div>
    );
  },
});
