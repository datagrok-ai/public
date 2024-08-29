import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { computed, defineComponent, PropType, ref, SlotsType, watchEffect } from 'vue';
import { IconFA } from '../elements';
import { useDraggable } from '@vueuse/core'

export const FoldableDialog = defineComponent({
  name: 'FoldableDialog',
  props: {
    title: {
      type: String,
      required: true,
    },
    foldable: {
      type: Object as PropType<true | false>,
      default: true
    }
  },
  slots: Object as SlotsType<{
    label?: any,
    default?: any,
  }>,
  emits: {
    closeClicked: () => {}
  },
  setup(props, {slots, emit}) {
    const isFolded = ref(false);

    const dialogEl = ref<HTMLElement | null>(null)

    // `style` will be a helper computed for `left: ?px; top: ?px;`
    const { x, y } = useDraggable(dialogEl, {
      initialValue: { x: 40, y: 40 },
    })

    const constrainedPosition = computed(() => {
      if (!dialogEl.value) return {x: x.value, y: y.value};

      const parentElement = dialogEl.value.parentElement ?? document.body;
      const xOffset = parentElement.getBoundingClientRect().left;
      const yOffset = parentElement.getBoundingClientRect().top;

      return {
        x: (xOffset < x.value) ? x.value: xOffset,
        y: (yOffset < y.value) ? y.value: yOffset,
      }
    });

    return () => (
      <div class='d4-dialog' 
        ref={dialogEl} 
        style={{
          'left': constrainedPosition.value.x + 'px',
          'top': constrainedPosition.value.y + 'px', 
          'position': 'fixed',
          'z-index': 3000,
        }}
      >
        <div class='flex items-center'>
          { props.foldable ? <IconFA 
            name={isFolded.value ? 'chevron-circle-right': 'chevron-circle-down'} 
            tooltip={isFolded.value ? 'Show contents': 'Hide contents' }
            onClick={() => isFolded.value = !isFolded.value}
            class='px-2'
          /> : null }
          { 
            <div class='d4-dialog-header justify-between gap-2'
              style={{'flex-grow': '1!important'}}
            >
              <div class='d4-dialog-title'> 
              { 
                slots.label ? 
                slots.label() : 
                props.title 
              }
              </div>
              <div class='flex items-end'>
                <IconFA
                  name='times'
                  tooltip='Close dialog'
                  onClick={() => emit('closeClicked')}
                />
              </div>
            </div> 
          }
        </div>
        <div class={{'d4-dialog-contents w-min': true, 'hidden': isFolded.value}}>
          { slots.default?.() }
        </div>
      </div>
    )
  }
});