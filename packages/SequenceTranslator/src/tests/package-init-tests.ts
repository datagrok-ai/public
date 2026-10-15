import {category, expectExceptionAsync, test} from '@datagrok-libraries/test/src/test';

import {OligoToolkitPackage} from '../apps/common/model/oligo-toolkit-package';

category('package-init', () => {
  test('seqHelperNotInitialized', async () => {
    const pkg = new OligoToolkitPackage();
    await expectExceptionAsync(async () => {pkg.seqHelper;},
      (e) => e instanceof Error && e.message.includes('.helmHelper is not initialized'));
  });
});
