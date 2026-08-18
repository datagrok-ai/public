import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {signal, computed, batch, Scope} from '@datagrok-libraries/u2';
import {buildDemo} from '../demo';

category('u2demo', () => {
  test('signals', async () => {
    const a = signal(1);
    const b = signal(2);
    const sum = computed(() => a.value + b.value);
    expect(sum.value, 3);
    batch(() => {
      a.value = 10;
      b.value = 32;
    });
    expect(sum.value, 42);
  });

  test('demo builds and disposes cleanly', async () => {
    const baseline = Scope.liveCount;
    const demo = buildDemo();
    expect(demo.root.childElementCount > 0, true);
    expect(Scope.liveCount > baseline, true);
    demo.dispose();
    expect(Scope.liveCount, baseline);
  });
});
