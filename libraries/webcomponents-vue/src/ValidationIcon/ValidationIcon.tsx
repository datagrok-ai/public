import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import type {ValidationIconT} from '@datagrok-libraries/webcomponents';
import type {ValidationIconInput} from '@datagrok-libraries/webcomponents/src/ValidationIcon/ValidationIcon';

declare global {
  namespace JSX {
    interface IntrinsicElements {
      'dg-validation-icon': ValidationIconT
    }
  }
}

export const ValidationIcon = Vue.defineComponent({
  name: 'ValidationIcon',
  props: {
    validationStatus: {
      type: Object as Vue.PropType<ValidationIconInput | undefined>,
      required: false
    },
    isScalar: {
      type: Boolean,
      default: true,
    }
  },
  emits: {
    consistencyReset: () => true,
    actionRequest: (_a: string) => true,
  },
  setup(props, {emit}) {
    const isScalar = Vue.computed(() => props.isScalar);
    const validationStatus = Vue.computed(() => props.validationStatus ? Vue.markRaw(props.validationStatus) : undefined);
    return () => (
      <dg-validation-icon
        isScalar={isScalar.value}
        validationStatus={validationStatus.value}
        onConsistencyReset={(_ev: CustomEvent) => emit('consistencyReset')}
        onActionRequest={(ev: CustomEvent<string>) => emit('actionRequest', ev.detail)}
      >
      </dg-validation-icon>
    );
  }
})
