import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

export const Overlapping = Vue.defineComponent({
  name: 'Overlapping',
  props: {
    isOverlapping: {
      type: Boolean,
      default: false,
    }
  },
  slots: Object as Vue.SlotsType<{
    default?: Vue.VNode[],
    overlapping?: Vue.VNode[],
  }>,
  setup(props, {slots}) {
    const styles = Vue.computed(() => {
      return props.isOverlapping ? {
        backgroundColor: '#dbdcdf',
        opacity: '0.8',
        position: 'absolute' as any,
        top: '0',
        right: '0',
        height: '100%',
        width: '100%',
        ...!slots.overlapping ? {
          alignItems: 'center',
          justifyContent: 'center',
          display: 'flex',
          flexDirection: 'column' as any
        }: {}
      } : {
        display: 'none',
      }
    })

    return () => <div>
      { slots.default?.() }
      <div style={styles.value}> 
        { slots.overlapping ? 
          slots.overlapping?.() : <Vue.Fragment>
            <label class={'ui-label'}> Updating... </label>
            <div class={'grok-loader'} style={{top: '10px', left: '-25px'}}> <div/> <div/> <div/> <div/> </div>
          </Vue.Fragment>
        }
      </div>
    </div>
  }
});