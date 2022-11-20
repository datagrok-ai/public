import {category, test, expect} from '@datagrok-libraries/utils/src/test';

category('ApiSamples', async () => {
  test('dummy', async () => {
    expect(1 == 1, true);
  })

  test('skipped dummy', async () => {
    expect(1 == 1, false);
  }, {skipReason: 'TASK-ID'})
})
