// Tests of demos parsing

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package-test';

import {category} from '@datagrok-libraries/test/src/test';
import {testTemplate} from './test-utils';
import {BIOREACTOR_MODEL_INFO} from '../demo/bioreactor';
import {PK_PD_MODEL_INFO} from '../demo/pk-pd';

// Demo models tests
category('Demo models', () => {
  testTemplate('bioreactor', BIOREACTOR_MODEL_INFO.equations);
  testTemplate('pk-pd', PK_PD_MODEL_INFO.equations);
}); // Demo models
