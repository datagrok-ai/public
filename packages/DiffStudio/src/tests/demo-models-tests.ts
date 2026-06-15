// Tests of demos parsing

import {category} from '@datagrok-libraries/test/src/test';
import {testIvpFile} from './test-utils';

// Demo models tests — equations come from the shipped `.ivp` files (single source of truth).
category('Demo models', () => {
  for (const name of ['acid-production', 'bioreactor', 'pk-pd', 'pollution'])
    testIvpFile(name, `${name}.ivp`);
}); // Demo models
