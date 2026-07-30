import * as Vue from 'vue';

const COL_GRIP = `url("data:image/svg+xml;utf8,<svg xmlns='http://www.w3.org/2000/svg' width='30' height='20'><path d='M 8 3 h 10 M 8 6 h 10 M 8 9 h 10' fill='none' stroke='%239497A0' stroke-width='1.25'/></svg>")`;
const ROW_GRIP = `url("data:image/svg+xml;utf8,<svg xmlns='http://www.w3.org/2000/svg' width='30' height='20'><path d='M 2 5 v 10 M 5 5 v 10 M 8 5 v 10' fill='none' stroke='%239497A0' stroke-width='1.25'/></svg>")`;

/**
 * Draggable divider that resizes an adjacent flex sibling, styled like ui.splitH/ui.splitV dividers.
 * Controlled component: bind the resized panel size via `size`/`update:size` (`v-model:size`).
 * `axis: 'x'` resizes width (col-resize), `axis: 'y'` resizes height (row-resize).
 * Set `reverse` when the resized panel is placed after the handle.
 */
export const ResizeHandle = Vue.defineComponent({
  name: 'ResizeHandle',
  props: {
    axis: {
      type: String as Vue.PropType<'x' | 'y'>,
      required: true,
    },
    size: {
      type: Number,
      required: true,
    },
    min: {
      type: Number,
      default: 0,
    },
    max: {
      type: Number,
      default: Number.POSITIVE_INFINITY,
    },
    reverse: {
      type: Boolean,
      default: false,
    },
  },
  emits: {
    'update:size': (_size: number) => true,
    'resizeEnd': (_size: number) => true,
  },
  setup(props, {emit}) {
    const onPointerDown = (e: PointerEvent) => {
      e.preventDefault();
      const handle = e.currentTarget as HTMLElement;
      handle.setPointerCapture(e.pointerId);
      const startPos = props.axis === 'x' ? e.clientX : e.clientY;
      const startSize = props.size;
      let lastSize = startSize;
      const onMove = (ev: PointerEvent) => {
        const delta = ((props.axis === 'x' ? ev.clientX : ev.clientY) - startPos) * (props.reverse ? -1 : 1);
        lastSize = Math.min(props.max, Math.max(props.min, startSize + delta));
        emit('update:size', lastSize);
      };
      const onUp = () => {
        handle.removeEventListener('pointermove', onMove);
        handle.releasePointerCapture(e.pointerId);
        emit('resizeEnd', lastSize);
      };
      handle.addEventListener('pointermove', onMove);
      handle.addEventListener('pointerup', onUp, {once: true});
    };

    return () => (
      <div
        style={{
          flexShrink: '0',
          background: `var(--grey-1) ${props.axis === 'x' ? COL_GRIP : ROW_GRIP} no-repeat center`,
          ...props.axis === 'x' ?
            {width: '5px', cursor: 'col-resize'} :
            {height: '5px', cursor: 'row-resize'},
        }}
        onPointerdown={onPointerDown}
      ></div>
    );
  },
});
