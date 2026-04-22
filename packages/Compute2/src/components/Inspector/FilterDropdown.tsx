import * as Vue from 'vue';

export interface FilterOption {
  value: string;
  label: string;
  detail?: string;
}

export const FilterDropdown = Vue.defineComponent({
  name: 'FilterDropdown',
  props: {
    options: {
      type: Array as Vue.PropType<FilterOption[]>,
      required: true,
    },
    modelValue: {
      type: Array as Vue.PropType<string[]>,
      required: true,
    },
    placeholder: {
      type: String,
      default: 'all',
    },
  },
  emits: {
    'update:modelValue': (_val: string[]) => true,
  },
  setup(props, {emit}) {
    const open = Vue.ref(false);
    const search = Vue.ref('');
    const rootRef = Vue.ref<HTMLElement>();

    const filtered = Vue.computed(() => {
      const q = search.value.toLowerCase();
      if (!q) return props.options;
      return props.options.filter((o) =>
        o.label.toLowerCase().includes(q) ||
        (o.detail?.toLowerCase().includes(q)));
    });

    const toggle = (value: string) => {
      const selected = props.modelValue;
      const idx = selected.indexOf(value);
      if (idx >= 0)
        emit('update:modelValue', selected.filter((v) => v !== value));
      else
        emit('update:modelValue', [...selected, value]);
    };

    const onClickOutside = (e: MouseEvent) => {
      if (rootRef.value && !rootRef.value.contains(e.target as Node))
        open.value = false;
    };
    Vue.onMounted(() => document.addEventListener('mousedown', onClickOutside));
    Vue.onUnmounted(() => document.removeEventListener('mousedown', onClickOutside));

    const triggerLabel = Vue.computed(() => {
      const n = props.modelValue.length;
      if (n === 0) return props.placeholder;
      if (n === 1) {
        const opt = props.options.find((o) => o.value === props.modelValue[0]);
        return opt?.label ?? props.modelValue[0];
      }
      return `${n} selected`;
    });

    return () => (
      <div ref={rootRef} style={{position: 'relative', display: 'inline-block'}}>
        <div
          style={{
            fontSize: '12px', padding: '1px 4px', cursor: 'pointer',
            border: '1px solid var(--border-color)', borderRadius: '2px',
            minWidth: '80px', userSelect: 'none', whiteSpace: 'nowrap',
          }}
          onClick={() => { open.value = !open.value; search.value = ''; }}
        >
          {triggerLabel.value}
          <span style={{marginLeft: '4px', fontSize: '10px'}}>&#9662;</span>
        </div>
        { open.value &&
          <div style={{
            position: 'absolute', zIndex: 1000, top: '100%', left: 0,
            border: '1px solid var(--border-color)', borderRadius: '2px',
            background: 'var(--white, #fff)', boxShadow: '0 2px 6px rgba(0,0,0,.15)',
            minWidth: '200px', maxHeight: '250px', display: 'flex',
            flexDirection: 'column',
          }}>
            <input
              type='text'
              placeholder='Filter...'
              v-model={search.value}
              style={{
                padding: '3px 6px', fontSize: '12px',
                border: 'none', borderBottom: '1px solid var(--border-color)', outline: 'none',
                background: 'transparent', color: 'inherit',
              }}
            />
            <div style={{overflow: 'auto', flex: 1}}>
              {filtered.value.map((opt) => (
                <label
                  key={opt.value}
                  style={{
                    display: 'flex', alignItems: 'center', gap: '4px',
                    padding: '2px 6px', fontSize: '12px', cursor: 'pointer',
                    whiteSpace: 'nowrap',
                  }}
                  onMouseenter={(e) => (e.currentTarget as HTMLElement).style.background = 'var(--steel-1, rgba(64, 96, 127, 0.1))'}
                  onMouseleave={(e) => (e.currentTarget as HTMLElement).style.background = ''}
                >
                  <input
                    type='checkbox'
                    checked={props.modelValue.includes(opt.value)}
                    onChange={() => toggle(opt.value)}
                    style={{flexShrink: 0}}
                  />
                  <span>{opt.label}</span>
                  { opt.detail &&
                    <span style={{opacity: 0.6, fontSize: '11px', marginLeft: '2px'}}>{opt.detail}</span>
                  }
                </label>
              ))}
              { filtered.value.length === 0 &&
                <div style={{padding: '4px 6px', fontSize: '11px', opacity: 0.6}}>No matches</div>
              }
            </div>
            { props.modelValue.length > 0 &&
              <div
                style={{
                  padding: '3px 6px', fontSize: '11px', cursor: 'pointer',
                  color: 'var(--blue-1)', borderTop: '1px solid var(--border-color)', textAlign: 'center',
                }}
                onClick={() => emit('update:modelValue', [])}
              >Clear all</div>
            }
          </div>
        }
      </div>
    );
  },
});
