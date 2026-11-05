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
    return () => (<button v-bind={attrs} onClick={(p) => emit('click', p)} is="dg-button">
      {slots.default ? slots.default() : ''}
    </button>);
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
    // Vue rewrites className whenever a class is bound, so the element's own classes are repeated here
    return () => (<button
      class={['ui-btn', 'ui-btn-ok', 'ui-btn-raised', {'d4-disabled': props.isDisabled}]}
      onClick={(p) => emit('click', p)}
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
      type: String as Vue.PropType<'spin' | 'pulse' | null>,
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
    const instance = Vue.getCurrentInstance();

    return () => {
      // a class binding, so a parent's class change keeps the disabled state
      return (<dg-icon-fa
        class={{'d4-disabled': props.isDisabled}}
        name={props.name}
        cursor={props.cursor}
        animation={props.animation}
        tooltip={props.tooltip}
        role={props.tooltip ? (instance?.vnode.props?.onClick ? 'button' : 'img') : undefined}
        aria-label={props.tooltip ?? undefined}
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
    tooltip: {
      type: String,
    },
  },
  emits: {
    'update:value': (_value: boolean) => true,
  },
  setup(props, {emit}) {
    return () => <dg-toggle-input
      caption={props.caption}
      value={props.value}
      tooltip={props.tooltip}
      onValueChanged={(event: any) => emit('update:value', event.detail)}
    />;
  },
});

export const ComboPopup = Vue.defineComponent({
  name: 'ComboPopup',
  props: {
    caption: {
      type: [String, HTMLElement] as Vue.PropType<string | HTMLElement>,
      required: true,
    },
    items: {
      type: Array as Vue.PropType<string[]>,
      default: () => [],
    },
    rightAligned: {
      type: Boolean,
      default: false,
    },
    tooltip: {
      type: String,
    },
  },
  emits: {
    'selected': (_value: {item: string, itemIdx: number}) => true,
  },
  setup(props, {emit}) {
    return () => <dg-combo-popup
      caption={props.caption}
      items={props.items}
      tooltip={props.tooltip}
      role={props.tooltip ? 'button' : undefined}
      aria-label={props.tooltip}
      rightAligned={props.rightAligned}
      onSelected={(e: any) => emit('selected', e.detail)}
    />;
  },
});
