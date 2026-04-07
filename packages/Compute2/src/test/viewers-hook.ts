import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {category, test, before, expect, awaitCheck} from '@datagrok-libraries/test/src/test';
import {useViewersHook} from '../composables/use-viewers-hook';
import {BehaviorSubject} from 'rxjs';
import * as Vue from 'vue';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';

category('Composables: Viewers hook', async () => {
  const mockStates = Vue.ref<Record<string, BehaviorSubject<Record<string, any>>> | undefined>();
  const mockCall = Vue.shallowRef(DG.Func.byName('Compute2:ObjectCooling2').prepare());

  before(async () => {
    mockStates.value = {
      ambTemp: Vue.markRaw(new BehaviorSubject<Record<string, any>>({})),
      simulation: Vue.markRaw(new BehaviorSubject<Record<string, any>>({})),
    };
    mockCall.value = DG.Func.byName('Compute2:ObjectCooling2').prepare();
  });

  test('Shoud run viewers hook on viewer changes', async () => {
    const mockData = grok.data.demo.demog();
    let lastIoName: string | undefined;
    let lastType: string | undefined;
    let lastViewer: DG.Viewer | undefined;
    let runCount = 0;
    const hook = (ioName: string, type: string, viewer?: DG.Viewer, meta?: any) => {
      lastIoName = ioName;
      lastType = type;
      lastViewer = viewer;
      runCount++;
    };
    const {setViewerRef} = useViewersHook(Vue.ref(hook), mockStates, mockCall as any);
    const grid = mockData.plot.grid();
    setViewerRef(grid, 'simulation', 'Grid');
    await awaitCheck(() => runCount > 0);
    expect(lastIoName, 'simulation');
    expect(lastType, 'Grid');
    expect(lastViewer, grid);

    const line = mockData.plot.line();
    setViewerRef(line, 'simulation', 'Line chart');
    await awaitCheck(() => runCount > 1);
    expect(lastIoName, 'simulation');
    expect(lastType, 'Line chart');
    expect(lastViewer, line);
  });

  test('Shoud run viewers hook on meta changes', async () => {
    const mockData = grok.data.demo.demog();
    let runCount = 0;
    let lastMeta: any = undefined;
    let lastIoName: string | undefined;
    let lastType: string | undefined;
    let lastViewer: DG.Viewer | undefined;
    const hook = (ioName: string, type: string, viewer?: DG.Viewer, meta?: any) => {
      // console.log(ioName, type, viewer, meta);
      lastIoName = ioName;
      lastType = type;
      lastViewer = viewer;
      lastMeta = meta;
      runCount++;
    };
    const {setViewerRef} = useViewersHook(Vue.ref(hook), mockStates, mockCall as any);
    const grid = mockData.plot.grid();
    setViewerRef(grid, 'simulation', 'Grid');
    await awaitCheck(() => runCount > 0);

    mockStates.value!.simulation.next({x: 1});
    await awaitCheck(() => runCount > 1);
    expect(lastIoName, 'simulation');
    expect(lastType, 'Grid');
    expect(lastViewer, grid);
    expectDeepEqual(lastMeta, {x: 1});
  });
});
