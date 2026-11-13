import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

import {tooltip} from './directives/tooltip';
import {onPress} from './on-press';

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
  emits: {
    click: (_e: MouseEvent) => true,
  },
  setup(_props, {slots, emit}) {
    // Vue rewrites className whenever a class is bound, so the element's own classes are repeated here
    return () => (<button class={['ui-btn', 'ui-btn-ok']} onClick={(p) => emit('click', p)} is="dg-button">
      {slots.default?.()}
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
    // why the button is disabled: its tooltip, and a warning when it is pressed anyway
    disabledReason: {
      type: String,
      required: false,
    },
  },
  emits: {
    click: (_e: MouseEvent) => true,
  },
  setup(props, {slots, emit}) {
    // a disabled button stays focusable, so it can tell why it does nothing; its tooltip shows only where the host
    // keeps pointer events on it, which d4-disabled turns off
    const onClick = (e: MouseEvent) => {
      if (!props.isDisabled)
        emit('click', e);
      else if (props.disabledReason)
        grok.shell.warning(props.disabledReason);
    };

    // Vue rewrites className whenever a class is bound, so the element's own classes are repeated here
    return () => Vue.withDirectives(<button
      class={['ui-btn', 'ui-btn-ok', 'ui-btn-raised', {'d4-disabled': props.isDisabled}]}
      aria-disabled={props.isDisabled ? 'true' : undefined}
      onClick={onClick}
      is="dg-big-button"
    >
      {slots.default?.()}
    </button>, [[tooltip, props.isDisabled ? props.disabledReason : undefined]]);
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
  emits: {
    click: (_e: MouseEvent) => true,
  },
  setup(props, {emit}) {
    const instance = Vue.getCurrentInstance();

    // a key press goes through a real click, so parents see the same event as from the mouse
    const pressOnKey = onPress((e) => (e.currentTarget as HTMLElement).click());

    return () => {
      const clickable = !!instance?.vnode.props?.onClick;
      const focusable = clickable && !props.isDisabled;
      // a class binding, so a parent's class change keeps the disabled state
      return (<dg-icon-fa
        class={{'d4-disabled': props.isDisabled}}
        name={props.name}
        cursor={props.cursor}
        animation={props.animation}
        tooltip={props.tooltip}
        role={props.tooltip ? (clickable ? 'button' : 'img') : undefined}
        aria-label={props.tooltip ?? undefined}
        tabindex={focusable ? 0 : undefined}
        faStyle={props.faStyle}
        onClick={(e: MouseEvent) => emit('click', e)}
        onKeydown={focusable ? pressOnKey : undefined}
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
