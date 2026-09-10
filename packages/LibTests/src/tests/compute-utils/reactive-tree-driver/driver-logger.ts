import {category, test} from '@datagrok-libraries/test/src/test';
import {DriverLogger} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/Logger';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';

category('ComputeUtils: Driver logger', async () => {
  test('Log is capped and keeps the newest entries', async () => {
    const logger = new DriverLogger();
    for (let i = 0; i < 5010; i++)
      logger.logLink('linkRunStarted', {linkUUID: `link-${i}`, prefix: [], id: `id-${i}`});
    const log = logger.logs$.value;
    expectDeepEqual(log.length, 5000, {prefix: 'Log length'});
    expectDeepEqual((log[log.length - 1] as any).linkUUID, 'link-5009', {prefix: 'Newest entry kept'});
    expectDeepEqual((log[0] as any).linkUUID, 'link-10', {prefix: 'Oldest entries dropped'});
  });
});
