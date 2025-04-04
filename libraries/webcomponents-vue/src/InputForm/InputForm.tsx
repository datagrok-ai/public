import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

import type {ValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';
import type {InputFormT} from '@datagrok-libraries/webcomponents';
import type {ConsistencyInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {injectInputBaseStatus, isInputInjected} from './utils';
import {BehaviorSubject, merge} from 'rxjs';
import {map, tap} from 'rxjs/operators';
import {useExtractedObservable} from '@vueuse/rxjs';

declare global {
  namespace JSX {
    interface IntrinsicElements {
      'dg-input-form': InputFormT
    }
  }
}

export const InputForm = Vue.defineComponent({
  name: 'InputForm',
  props: {
    funcCall: {
      type: Object as Vue.PropType<DG.FuncCall>,
      required: true,
    },
    skipInit: {
      type: Boolean,
      default: true,
    },
    validationStates: {
      type: Object as Vue.PropType<Record<string, ValidationResult>>,
    },
    consistencyStates: {
      type: Object as Vue.PropType<Record<string, ConsistencyInfo>>,
    },
    callMeta: {
      type: Object as Vue.PropType<Record<string, BehaviorSubject<any>>>,
    },
    isReadonly: {
      type: Boolean,
    },
  },
  emits: {
    formReplaced: (_a: DG.InputForm | undefined) => true,
    inputChanged: (_a: DG.EventData<DG.InputArgs>) => true,
    validationChanged: (_a: boolean) => true,
    actionRequested: (_actionUuid: string) => true,
    consistencyReset: (_ioName: string) => true,
  },
  setup(props, {emit}) {
    const currentCall = Vue.computed(() => Vue.markRaw(props.funcCall));
    const validationStates = Vue.computed(() => props.validationStates);
    const consistencyStates = Vue.computed(() => props.consistencyStates);
    const isReadonly = Vue.computed(() => props.isReadonly);
    const skipInit = Vue.computed(() => props.skipInit);

    const states = Vue.reactive({
      meta: {} as Record<string, any>,
    });

    useExtractedObservable(() => props.callMeta, (meta) => {
      states.meta = {};
      const entries = Object.entries(meta).map(([name, state$]) => state$.pipe(map((s) => [name, s] as const)));
      return merge(...entries).pipe(
        tap(([k, val]) => {
          states.meta[k] = val ? Vue.markRaw(val) : undefined;
        }),
      );
    });

    const currentForm = Vue.shallowRef(undefined as undefined | DG.InputForm);

    Vue.watch([currentCall, currentForm, validationStates, consistencyStates], ([call, form, validationStates, consistencyStates]) => {
      if (!form || !call) return;

      [...call.inputParams.values()]
        .filter((param) => (form.getInput(param.property.name)))
        .forEach((param) => {
          const input = form.getInput(param.property.name);
          const validation = validationStates?.[param.property.name];
          const consistency = consistencyStates?.[param.property.name];

          if (!isInputInjected(input))
            injectInputBaseStatus(emit, param.property.name, input);

          (input as any).setStatus({validation, consistency});

          input.enabled = !isReadonly.value && consistency?.restriction !== 'disabled';
        });
    });

    Vue.watch([currentCall, currentForm, states.meta], ([call, form, meta]) => {
      if (!form || !call) return;

      [...call.inputParams.values()]
        .filter((param) => (form.getInput(param.property.name)))
        .forEach((param) => {
          const input = form.getInput(param.property.name);
          const paramItems = meta[param.property.name]?.['items'];
          if (input.inputType === DG.InputType.Choice) {
            input.notify = false;
            const currentValue = param.value;
            try {
              if (paramItems)
                (input as DG.ChoiceInput<any>).items = paramItems;
              else if (param.property.options.choices)
                (input as DG.ChoiceInput<any>).items = JSON.parse(param.property.options.choices);
            } catch(e) {
              console.error(e);
            } finally {
              input.value = currentValue;
              input.notify = true;
            }
          }
          // TODO: FormApi
          // const rangeMeta = meta[param.property.name]?.['range'];
          //if (rangeMeta && (input.inputType === DG.InputType.Float || input.inputType === DG.InputType.Int)) {}
          const hideMeta = meta[param.property.name]?.['hidden'];
          if (hideMeta)
            input.root.style.display = 'none';
          else
            input.root.style.display = 'flex';
        });
    });

    const formReplacedCb = async (event: {detail?: DG.InputForm}) => {
      const form = event.detail ? Vue.markRaw(event.detail) : undefined
      currentForm.value = form;
      emit('formReplaced', form);
    };

    return () =>
      <dg-input-form
        skipInit={skipInit.value}
        funcCall={currentCall.value}
        onFormReplaced={formReplacedCb}
        onInputChanged={(ev: CustomEvent<DG.EventData<DG.InputArgs>>) => emit('inputChanged', ev.detail)}
        onValidationChanged={(ev: CustomEvent<boolean>) => emit('validationChanged', ev.detail)}>
      </dg-input-form>;
  },
});
