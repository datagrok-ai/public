import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

import {useObservable, useSubscription} from '@vueuse/rxjs'

export const DockedPanel = Vue.defineComponent({
  name: 'DockedPanel',
  props: {
    dockType: {
      type: String as Vue.PropType<DG.DockType>,
    },
    title: {
      type: String,
    },
    ratio: {
      type: Number
    }
  },
  slots: Object as Vue.SlotsType<{
    default?: any,
  }>,
  emits: {
    closed: () => {}
  },
  setup(props, {slots, emit}) {
    const element = Vue.ref(null as null | HTMLElement)
    let parentView = null as null | DG.ViewBase;
    let dockNode = null as null | DG.DockNode;

    Vue.onMounted(async () => {
      await Vue.nextTick();
      parentView = grok.shell.v;

      if (!element.value) return;

      dockNode = grok.shell.dockManager.dock(element.value, props.dockType, null, props.title, props.ratio) 
    })

    const closedPanel = useObservable(grok.shell.dockManager.onClosed);
    Vue.watch(closedPanel, (val) => {
      if (val === element.value) {
        emit('closed')
      }
    })

    const changedView = useObservable(grok.events.onCurrentViewChanged);
    Vue.watch(changedView, (newView) => {
      if (newView != parentView) {
        if (!dockNode) return;

        grok.shell.dockManager.close(dockNode);
        emit('closed')
      } else {
        if (!element.value) return;

        dockNode = grok.shell.dockManager.dock(element.value, props.dockType, null, props.title, props.ratio) 
      }
    })

    Vue.onUnmounted(() => {
      if (!dockNode) return;

      grok.shell.dockManager.close(dockNode);
    });

    return () => 
      <div ref={element} style={{width: '100%', height: '100%'}}> { slots.default?.() } </div>
  }
});