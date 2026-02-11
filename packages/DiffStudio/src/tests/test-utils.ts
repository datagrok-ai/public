// Tests of equations parsing

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package-test';

import {expect, test} from '@datagrok-libraries/test/src/test';

import {IVP, getIVP, getScriptLines, getScriptParams} from '../scripting-tools';
import {DF_NAME} from '../constants';

const TIMEOUT = 10000;
const MIN_ROWS = 1;

/** Template for testing solving ODEs */
export function testTemplate(name: string, problem: string): void {
  test(name, async () => {
    let ivp: IVP | null = null;
    let code: string | null = null;
    let script: DG.Script | null = null;
    let params: Record<string, number> | null = null;
    let call: DG.FuncCall | null = null;
    let solution: DG.DataFrame | null = null;

    console.log(problem);

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
