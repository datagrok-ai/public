import * as Vue from 'vue';

const COL_GRIP = `url("data:image/svg+xml;utf8,<svg xmlns='http://www.w3.org/2000/svg' width='30' height='20'><path d='M 8 3 h 10 M 8 6 h 10 M 8 9 h 10' fill='none' stroke='%239497A0' stroke-width='1.25'/></svg>")`;
const ROW_GRIP = `url("data:image/svg+xml;utf8,<svg xmlns='http://www.w3.org/2000/svg' width='30' height='20'><path d='M 2 5 v 10 M 5 5 v 10 M 8 5 v 10' fill='none' stroke='%239497A0' stroke-width='1.25'/></svg>")`;

interface GhostState {
  pos: number;
  crossStart: number;
  crossSize: number;
}

/**
 * Draggable divider that resizes an adjacent flex sibling, styled like ui.splitH/ui.splitV dividers.
 * Controlled component: bind the resized panel size via `size`/`update:size` (`v-model:size`).
 * `axis: 'x'` resizes width (col-resize), `axis: 'y'` resizes height (row-resize).
 * Set `reverse` when the resized panel is placed after the handle.
 * By default drags move a ghost line and `update:size` fires once on release, since resizing
 * canvas viewers per frame is expensive; set `live` for cheap panels that can follow the drag.
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
    live: {
      type: Boolean,
      default: false,
    },
  },
  emits: {
    'update:size': (_size: number) => true,
    'resizeEnd': (_size: number) => true,
  },
  setup(props, {emit}) {
    const ghost = Vue.ref<GhostState | null>(null);

    const onPointerDown = (e: PointerEvent) => {
      e.preventDefault();
      const handle = e.currentTarget as HTMLElement;
      handle.setPointerCapture(e.pointerId);
      const rect = handle.getBoundingClientRect();
      const handleStart = props.axis === 'x' ? rect.left : rect.top;
      const crossStart = props.axis === 'x' ? rect.top : rect.left;
      const crossSize = props.axis === 'x' ? rect.height : rect.width;
      const startPos = props.axis === 'x' ? e.clientX : e.clientY;
      const startSize = props.size;
      const direction = props.reverse ? -1 : 1;
      let lastSize = startSize;
      let pendingFrame = 0;

      if (!props.live)
        ghost.value = {pos: handleStart, crossStart, crossSize};

      const onMove = (ev: PointerEvent) => {
        const delta = ((props.axis === 'x' ? ev.clientX : ev.clientY) - startPos) * direction;
        lastSize = Math.min(props.max, Math.max(props.min, startSize + delta));
        if (props.live) {
          // coalesce to one emit per animation frame
          if (!pendingFrame) {
            pendingFrame = requestAnimationFrame(() => {
              pendingFrame = 0;
              emit('update:size', lastSize);
            });
          }
        } else {
          ghost.value = {pos: handleStart + (lastSize - startSize) * direction, crossStart, crossSize};
        }
      };
      const onUp = () => {
        handle.removeEventListener('pointermove', onMove);
        handle.releasePointerCapture(e.pointerId);
        if (pendingFrame) {
          cancelAnimationFrame(pendingFrame);
          pendingFrame = 0;
        }
        ghost.value = null;
        emit('update:size', lastSize);
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
      >
        { ghost.value != null &&
          <div style={{
            position: 'fixed', zIndex: '1000', pointerEvents: 'none',
            background: 'var(--blue-1, #2083d5)', opacity: '0.5',
            ...props.axis === 'x' ?
              {
                left: `${ghost.value.pos}px`, width: '3px',
                top: `${ghost.value.crossStart}px`, height: `${ghost.value.crossSize}px`,
              } :
              {
                top: `${ghost.value.pos}px`, height: '3px',
                left: `${ghost.value.crossStart}px`, width: `${ghost.value.crossSize}px`,
              },
          }}></div>
        }
      </div>
    );
  },
});
