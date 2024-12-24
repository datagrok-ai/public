import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

export const FunctionCard = Vue.defineComponent({
  name: 'FunctionCard',
  props: {
    func: {
      type: Object as Vue.PropType<DG.Func>,
      required: true,
    },
  },
  emits: {
    'click': () => {},
  },
  setup(props, {emit}) {
    const currentFunc = Vue.computed(() => props.func);

    const language = Vue.computed(() =>
      currentFunc.value instanceof DG.Script ?
        currentFunc.value.language: 'javascript',
    );

    const iconBackground = Vue.computed(() => `background-image: url("/images/entities/${language.value}.png");`);

    return () => <div
      onClick={() => emit('click')}
      class={'grok-gallery-grid-item-wrapper'}
      style={{cursor: 'pointer'}}
    >
      <div class={'grok-gallery-grid-item grok-scripting-script d4-flex-col d4-gallery-card entity-script'}>
        <div class={'d4-flex-col'}>
          <span class={'d4-link-label'}>
            <i style={iconBackground.value} class={'grok-icon image-icon'}></i>
            <label class={'grok-gallery-grid-item-title'}>
              {currentFunc.value.friendlyName}
            </label>
          </span>
          <label id="description">
            {currentFunc.value.description.length > 0 ? currentFunc.value.description: 'No description provided'}
          </label>
        </div>
      </div>
    </div>;
  },
});
