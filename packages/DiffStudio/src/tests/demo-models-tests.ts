// Tests of demos parsing

import {category} from '@datagrok-libraries/test/src/test';
import {testIvpFile} from './test-utils';

// Demo models tests — equations come from the shipped `.ivp` files (single source of truth).
category('Demo models', () => {
  for (const name of ['acid-production', 'bioreactor', 'pk-pd', 'pollution'])
    testIvpFile(name, `library/${name}.ivp`);

  // Model-Hub / demo models — divergent from the library examples (own equations).
  testIvpFile('pk-pd model', 'models/pk-pd.ivp');
  testIvpFile('bioreactor model', 'models/bioreactor.ivp');
}); // Demo models
