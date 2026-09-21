import {category, test, before} from '@datagrok-libraries/test/src/test';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {TestScheduler} from 'rxjs/testing';
import {of} from 'rxjs';
import {delay, mapTo} from 'rxjs/operators';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {createTestScheduler} from '../../../test-utils';

const mixedSpeedHandlerConfig: PipelineConfiguration = {
  id: 'root',
  type: 'static',
  steps: [
    {id: 'step1', nqName: 'LibTests:TestAdd2'},
    {id: 'step2', nqName: 'LibTests:TestMul2'},
  ],
  links: [{
    id: 'mixed-link',
    from: 'in1:step1/b',
    to: 'out1:step2/a',
    handler({controller}) {
      const v = controller.getFirst('in1');
      controller.setAll('out1', v);
      if (v === 'slow')
        return of(null).pipe(delay(100), mapTo(void(0)));
    },
  }],
};

category('ComputeUtils: Driver link running state', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Run scheduled right after a finished run reports running', async () => {
    const pconf = await getProcessedConfig(mixedSpeedHandlerConfig);
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const store = tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore();
      const link = [...tree.linksState.links.values()].find((l) => l.matchInfo.spec.id === 'mixed-link')!;
      let isRunning: boolean | undefined;
      link.isRunning$.subscribe((val) => isRunning = val);
      cold('-a').subscribe(() => {
        // a sync run finishes and an async run is scheduled within the same frame,
        // so wall-clock/frame stamps of both are equal
        store.setState('b', 'fast');
        link.trigger();
        store.setState('b', 'slow');
        link.trigger();
      });
      cold('--------a').subscribe(() => {
        expectDeepEqual(isRunning, true, {prefix: 'Running while the async handler is in flight'});
      });
      cold('150ms a').subscribe(() => {
        expectDeepEqual(isRunning, false, {prefix: 'Not running after the async handler finished'});
      });
    });
  });

  test('Async run reports running for its whole duration', async () => {
    const pconf = await getProcessedConfig(mixedSpeedHandlerConfig);
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const store = tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore();
      const link = [...tree.linksState.links.values()].find((l) => l.matchInfo.spec.id === 'mixed-link')!;
      let isRunning: boolean | undefined;
      link.isRunning$.subscribe((val) => isRunning = val);
      cold('-a').subscribe(() => {
        store.setState('b', 'slow');
      });
      cold('50ms a').subscribe(() => {
        expectDeepEqual(isRunning, true, {prefix: 'Running mid-flight'});
      });
      cold('150ms a').subscribe(() => {
        expectDeepEqual(isRunning, false, {prefix: 'Not running after completion'});
      });
    });
  });
});
