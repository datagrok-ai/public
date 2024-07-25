import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {defineComponent, ref, watch} from 'vue';
import type {DGBigButtonT, DGButtonT, DGIconFAT, DGSplitH} from '@datagrok-libraries/webcomponents/src';

declare global {
  namespace JSX {
    interface IntrinsicElements {
      'dg-button': DGButtonT,
      'dg-big-button': DGBigButtonT,
      'dg-split-h': DGSplitH,
      'dg-icon-fa': DGIconFAT,
    }
  }
}

export const Button = defineComponent({
  name: 'Button',
  emits: ['click'],
  setup(_props, {slots, attrs, emit}) {
    return () => (<button v-bind={attrs} onClick={(p) => emit('click', p)} is="dg-button">{slots.default ? slots.default() : ''}</button>);
  },
});

export const BigButton = defineComponent({
  name: 'BigButton',
  emits: ['click'],
  setup(_props, {slots, attrs, emit}) {
    return () => (<button v-bind={attrs} onClick={(p) => emit('click', p)} is="dg-big-button">{slots.default ? slots.default() : ''}</button>);
  },
});


export const SplitH = defineComponent({
  name: 'SplitH',
  props: {
    resize: Boolean,
  },
  setup(props, {slots, attrs, emit}) {
    return () =>{
      return (<dg-split-h
        resize={props.resize}
        v-bind={attrs}
        style={{height: '100%'}}
      >
        {slots.default ? slots.default() : []}
      </dg-split-h>);
    };
  },
});

export const IconFA = defineComponent({
  name: 'IconFA',
  props: {
    name: String,
    cursor: {
      type: String,
      default: 'pointer',
    },
  },
  emits: [
    'click',
  ],
  setup(props, {emit}) {
    return () => {
      return (<dg-icon-fa
        name={props.name}
        cursor={props.cursor}
        onClick={(e: Event) => emit('click', e)}
      >
      </dg-icon-fa>);
    };
  },
});
