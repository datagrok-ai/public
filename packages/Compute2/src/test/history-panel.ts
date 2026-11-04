import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {provideDgViewService} from '@datagrok-libraries/webcomponents-vue';
import {useHistoryPanel} from '../composables/use-history-panel';

type Panel = ReturnType<typeof useHistoryPanel>;
// position 0 and 1 stand for the RFV and the PipelineView, siblings in TreeWizard
type Step = {uuid: string, nqName: string, position: 0 | 1};

const StepView = Vue.defineComponent({
  props: {nqName: {type: String, required: true}},
  emits: {ready: (_panel: Panel) => true},
  setup(props, {emit}) {
    emit('ready', useHistoryPanel(() => props.nqName));
    return () => Vue.h('div');
  },
});

// keyed by uuid like TreeWizard, so a load or a navigation rebuilds the step view
function mountSteps(view: DG.ViewBase, first: Step) {
  const shown = Vue.shallowRef(first);
  let panel: Panel | undefined;
  const app = Vue.createApp({
    setup: () => () => Vue.h('div', ([0, 1] as const).map((position) => shown.value.position === position ?
      Vue.h(StepView, {key: shown.value.uuid, nqName: shown.value.nqName, onReady: (ready: Panel) => panel = ready}) :
      null)),
  });
  provideDgViewService(app, view);
  app.mount(ui.div());
  return {
    panel: () => panel!,
    show: async (step: Step) => {
      shown.value = step;
      await Vue.nextTick();
    },
    unmount: () => app.unmount(),
  };
}

const open = (panel: Panel, readonly: boolean) => {
  panel.historyReadonly.value = readonly;
  panel.historyHidden.value = false;
};

category('History panel', () => {
  test('Stays open when a load rebuilds the step view', async () => {
    const steps = mountSteps(DG.View.create(), {uuid: '1', nqName: 'Pkg:Step', position: 0});
    try {
      open(steps.panel(), true);
      await steps.show({uuid: '2', nqName: 'Pkg:Step', position: 0});
      expect(steps.panel().historyHidden.value, false);
      expect(steps.panel().historyReadonly.value, true);

      await steps.show({uuid: '3', nqName: 'Pkg:Step', position: 0});
      expect(steps.panel().historyHidden.value, false);
    } finally {
      steps.unmount();
    }
  });

  test('Stays closed when a load rebuilds the step view', async () => {
    const steps = mountSteps(DG.View.create(), {uuid: '1', nqName: 'Pkg:Step', position: 0});
    try {
      await steps.show({uuid: '2', nqName: 'Pkg:Step', position: 0});
      expect(steps.panel().historyHidden.value, true);
    } finally {
      steps.unmount();
    }
  });

  test('Closes on a step of another function and stays closed back', async () => {
    const steps = mountSteps(DG.View.create(), {uuid: '1', nqName: 'Pkg:First', position: 0});
    try {
      open(steps.panel(), false);
      await steps.show({uuid: '2', nqName: 'Pkg:Second', position: 0});
      expect(steps.panel().historyHidden.value, true);

      await steps.show({uuid: '1', nqName: 'Pkg:First', position: 0});
      expect(steps.panel().historyHidden.value, true);
    } finally {
      steps.unmount();
    }
  });

  test('Closes between a step and a pipeline both ways', async () => {
    const steps = mountSteps(DG.View.create(), {uuid: '1', nqName: 'Pkg:Pipeline', position: 1});
    try {
      open(steps.panel(), false);
      await steps.show({uuid: '2', nqName: 'Pkg:Step', position: 0});
      expect(steps.panel().historyHidden.value, true);
      await steps.show({uuid: '1', nqName: 'Pkg:Pipeline', position: 1});
      expect(steps.panel().historyHidden.value, true);

      open(steps.panel(), false);
      await steps.show({uuid: '2', nqName: 'Pkg:Step', position: 0});
      open(steps.panel(), false);
      await steps.show({uuid: '1', nqName: 'Pkg:Pipeline', position: 1});
      expect(steps.panel().historyHidden.value, true);
      await steps.show({uuid: '2', nqName: 'Pkg:Step', position: 0});
      expect(steps.panel().historyHidden.value, true);
    } finally {
      steps.unmount();
    }
  });

  test('Is not handed to another view', async () => {
    const first = mountSteps(DG.View.create(), {uuid: '1', nqName: 'Pkg:Step', position: 0});
    open(first.panel(), false);
    first.unmount();
    const second = mountSteps(DG.View.create(), {uuid: '2', nqName: 'Pkg:Step', position: 0});
    try {
      expect(second.panel().historyHidden.value, true);
    } finally {
      second.unmount();
    }
  });
});
