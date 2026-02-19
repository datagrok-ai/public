// Tests of equations parsing

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package-test';

import {category} from '@datagrok-libraries/test/src/test';

import {TEMPLATES, ENERGY_N_CONTROL} from '../templates';
import {USE_CASES} from '../use-cases';
import {testTemplate} from './test-utils';

// Correctness tests
category('Features', () => {
  testTemplate('Basic project', TEMPLATES.BASIC);
  testTemplate('Project structs', TEMPLATES.ADVANCED);
  testTemplate('Annotating params', TEMPLATES.EXTENDED);
  testTemplate('Cyclic process', USE_CASES.PK_PD);
  testTemplate('Multistage model', USE_CASES.ACID_PROD);
  testTemplate('Output control', USE_CASES.NIMOTUZUMAB);
  testTemplate('Value lookups', USE_CASES.NIMOTUZUMAB);
  testTemplate('Output expressions & use of JS code in model', ENERGY_N_CONTROL);
}); // Correctness
