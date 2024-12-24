import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

import type {DGBigButtonT, DGButtonT, DGComboPopupT, DGIconFAT, DGToggleInputT} from '@datagrok-libraries/webcomponents';

declare global {
  namespace JSX {
    interface IntrinsicElements {
      'dg-button': DGButtonT,
      'dg-big-button': DGBigButtonT,
      'dg-icon-fa': DGIconFAT,
      'dg-toggle-input': DGToggleInputT
      'dg-combo-popup': DGComboPopupT
    }
  }
}

export const Button = Vue.defineComponent({
  name: 'Button',
  emits: ['click'],
  setup(_props, {slots, attrs, emit}) {
    return () => (<button v-bind={attrs} onClick={(p) => emit('click', p)} is="dg-button">{slots.default ? slots.default() : ''}</button>);
  },
});

export const BigButton = Vue.defineComponent({
  name: 'BigButton',
  props: {
    isDisabled: {
      type: Boolean,
      default: false,
    },
  },
  emits: ['click'],
  setup(props, {slots, emit}) {
    const btn = Vue.shallowRef<HTMLElement | null>(null);

    Vue.watch([() => props.isDisabled, btn] as const, ([isDisabled, btn]) => {
      if (!btn) return;

      btn.classList.toggle('d4-disabled', isDisabled);
    }, {immediate: true});

    return () => (<button
      onClick={(p) => emit('click', p)}
      ref={btn}
      is="dg-big-button"
    >
      {slots.default ? slots.default() : ''}
    </button>);
  },
});

export const IconFA = Vue.defineComponent({
  name: 'IconFA',
  props: {
    name: String,
    cursor: {
      type: String,
      default: 'pointer',
    },
    animation: {
      type: Object as Vue.PropType<'spin' | 'pulse' | null>,
      default: null,
    },
    tooltip: {
      type: String as Vue.PropType<string | null>,
      default: null,
    },
    faStyle: {
      type: String as Vue.PropType<'fal' | 'fas' | 'far' | 'fad'>,
      default: 'fal',
    },
    isDisabled: {
      type: Boolean,
      default: false,
    },
  },
  emits: [
    'click',
  ],
  setup(props, {emit}) {
    const icon = Vue.shallowRef<HTMLElement | null>(null);

    Vue.watch([() => props.isDisabled, icon] as const, ([isDisabled, btn]) => {
      if (!btn) return;

      btn.classList.toggle('d4-disabled', isDisabled);
    }, {immediate: true});

    return () => {
      return (<dg-icon-fa
        name={props.name}
        ref={icon}
        cursor={props.cursor}
        animation={props.animation}
        tooltip={props.tooltip}
        faStyle={props.faStyle}
        onClick={(e: Event) => emit('click', e)}
      >
      </dg-icon-fa>);
    };
  },
});

export const ToggleInput = Vue.defineComponent({
  name: 'ToggleInput',
  props: {
    caption: {
      type: String,
      required: true,
    },
    value: {
      type: Boolean,
      default: false,
    },
    nullable: {
      type: Boolean,
      default: false,
    },
    tooltip: {
      type: String,
    },
  },
  emits: {
    'update:value': (value: boolean) => value,
  },
  setup(props, {emit, attrs}) {
    return () => <dg-toggle-input
      caption={props.caption}
      value={props.value}
      nullable={props.nullable}
      tooltip={props.tooltip}
      onValueChanged={(event: any) => emit('update:value', event.detail)}
    />;
  },
});

export const ComboPopup = Vue.defineComponent({
  name: 'ComboPopup',
  props: {
    caption: {
      type: Object as Vue.PropType<String | HTMLElement>,
      required: true,
    },
    items: {
      type: Array<String>,
      default: [],
    },
    tooltip: {
      type: String,
    },
  },
  emits: {
    'selected': (value: {item: string, itemIdx: number}) => value,
  },
  setup(props, {emit}) {
    return () => <dg-combo-popup
      //@ts-ignore
      caption={props.caption}
      items={props.items}
      tooltip={props.tooltip}
      onSelected={(e: any) => emit('selected', e.detail)}
    />;
  },
});
