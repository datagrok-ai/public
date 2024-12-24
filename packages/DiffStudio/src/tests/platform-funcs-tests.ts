// Tests of platform functions

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package-test';

import {category, expect, test} from '@datagrok-libraries/utils/src/test';

import {USE_CASES} from '../use-cases';
import {IVP} from '../scripting-tools';

// Platform functions tests
category('Functions', () => {
  test('Model-to-script', async () => {
    let ivp: IVP | null = null;
    let msg: string;

    try {
      const serializer = await grok.functions.eval('DiffStudio:serializeEquations');
      const call = serializer.prepare({problem: USE_CASES.BIOREACTOR});
      await call.call();
      ivp = call.getParamValue('serialization');
    } catch (err) {
      msg = (err instanceof Error) ? err.message : 'platform issue';
    }

    expect((ivp !== null), true, `Serialization failed: ${msg}`);

    let scriptCode: string | null = null;

    try {
      const coder = await grok.functions.eval('DiffStudio:odesToCode');
      const call = coder.prepare({serialization: ivp});
      await call.call();
      scriptCode = call.getParamValue('code');
    } catch (err) {
      msg = (err instanceof Error) ? err.message : 'platform issue';
    }

    expect((scriptCode !== null), true, `Generating JS-code failed: ${msg}`);
  });

  test('Solve equations', async () => {
    let solution: DG.DataFrame | null = null;
    let msg: string | null = null;

    try {
      const solver = await grok.functions.eval('DiffStudio:solveODE');
      const call = solver.prepare({problem: USE_CASES.ACID_PROD});
      await call.call();
      solution = call.getParamValue('solution');
    } catch (err) {
      msg = (err instanceof Error) ? err.message : 'platform issue';
    }

    expect((solution !== null), true, `Solving equations failed: ${msg}`);
  });
}); // Correctness
