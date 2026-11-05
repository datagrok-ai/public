import * as ui from 'datagrok-api/ui';
import * as Vue from 'vue';

const containerByApp = new WeakMap<Vue.App, HTMLElement>();

/** Like KeepAlive, but the kept instance outlives the parent: a rebuilt parent rendering content with the same key
 * gets it back. One instance per Vue app; it is rendered detached and moved into the current placeholder. */
export const KeptAlive = Vue.defineComponent({
  name: 'KeptAlive',
  setup(_props, {slots}) {
    const appContext = Vue.getCurrentInstance()!.appContext;
    const app = appContext.app;
    const placeholder = Vue.shallowRef<HTMLElement | undefined>(undefined);
    let content: Vue.VNode | undefined;

    const containerForApp = () => {
      let container = containerByApp.get(app);
      if (!container) {
        container = ui.div([], 'h-full');
        containerByApp.set(app, container);
        app.onUnmount(() => Vue.render(null, container!));
      }
      return container;
    };

    const attach = () => {
      const container = containerForApp();
      content!.appContext = appContext;
      Vue.render(content!, container);
      if (container.parentElement !== placeholder.value)
        placeholder.value!.append(container);
    };
    Vue.onMounted(attach);
    Vue.onUpdated(attach);
    Vue.onBeforeUnmount(() => {
      const container = containerByApp.get(app);
      if (container?.parentElement === placeholder.value)
        container!.remove();
    });

    return () => {
      content = slots.default!()[0];
      return Vue.h('div', {ref: placeholder});
    };
  },
});
