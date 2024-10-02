// Tests of equations parsing

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package-test';

import {category, expect, test, timeout} from '@datagrok-libraries/utils/src/test';

import {TEMPLATES, ENERGY_N_CONTROL} from '../templates';
import {USE_CASES} from '../use-cases';
import {IVP, getIVP, getScriptLines, getScriptParams} from '../scripting-tools';
import {DF_NAME} from '../constants';

const TIMEOUT = 10000;
const MIN_ROWS = 1;

/** Template for testing solving ODEs */
const testTeamplate = (name: string, problem: string) => {
  test(name, async () => {
    let ivp: IVP | null = null;
    let code: string | null = null;
    let script: DG.Script | null = null;
    let params: Record<string, number> | null = null;
    let call: DG.FuncCall | null = null;
    let solution: DG.DataFrame | null = null;

    try {
      ivp = getIVP(problem);
      code = getScriptLines(ivp).join('\n');
      script = DG.Script.create(code);
      params = getScriptParams(ivp);
      call = script.prepare(params);

      await call.call();

      solution = call.outputs[DF_NAME];
    } catch (e) {}

    expect( ivp !== null, true, 'Failed to parse equations');
    expect( code !== null, true, 'Failed to create JS-code');
    expect( script !== null, true, 'Failed to create the platform script');
    expect( params !== null, true, 'Failed to get script params');
    expect( call !== null, true, 'Failed to get funccall');
    expect( solution !== null, true, 'Failed to apply numerical method');
    expect( solution.rowCount >= MIN_ROWS, true, 'Solving failed: an empty dataframe');
  }, {timeout: TIMEOUT});
};

// Correctness tests
category('Features', () => {
  testTeamplate('Basic project', TEMPLATES.BASIC);
  testTeamplate('Project structs', TEMPLATES.ADVANCED);
  testTeamplate('Annotating params', TEMPLATES.EXTENDED);
  testTeamplate('Cyclic process', USE_CASES.PK_PD);
  testTeamplate('Multistage model', USE_CASES.ACID_PROD);
  testTeamplate('Output control', USE_CASES.NIMOTUZUMAB);
  testTeamplate('Value lookups', USE_CASES.NIMOTUZUMAB);
  testTeamplate('Output expressions & use of JS code in model', ENERGY_N_CONTROL);
}); // Correctness
