import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {Scope} from '@datagrok-libraries/u2';
import {funcForm} from '@datagrok-libraries/u2/dg';

const CAP = 500;
const SETTLE_MS = 5000;

/** The functions with inputs, spread over every package when there are more than {@link CAP}. */
function sweptFunctions(): DG.Func[] {
  const all = DG.Func.find().filter((f) => f.inputs.length > 0).sort((a, b) => a.nqName.localeCompare(b.nqName));
  const step = Math.ceil(all.length / CAP);
  return all.filter((_, i) => i % step === 0);
}

const text = (args: unknown[]) => args.map((a) => a instanceof Error ? a.message : String(a)).join(' ');

// B-3: one open table that lost its name broke every function's form
category('u2/dg: funcForm sweep', () => {
  test('every function with inputs builds its form, no throw and no console error', async () => {
    const funcs = sweptFunctions();
    const named = grok.shell.addTable(grok.data.demo.demog(20));
    const nameless = grok.shell.addTable(DG.DataFrame.fromCsv('a\n1'));
    nameless.name = null as unknown as string;
    const thrown: string[] = [];
    const errors: string[] = [];
    const problems: string[] = [];
    let current = '';
    let unprepared = 0;
    let unsettled = 0;
    const consoleError = console.error;
    console.error = (...args: unknown[]) => {
      errors.push(`${current}: ${text(args)}`);
      consoleError.apply(console, args);
    };
    const onError = (e: ErrorEvent) => errors.push(`${current}: ${e.message}`);
    const onRejection = (e: PromiseRejectionEvent) => errors.push(`${current}: ${text([e.reason])}`);
    window.addEventListener('error', onError);
    window.addEventListener('unhandledrejection', onRejection);
    try {
      for (const func of funcs) {
        current = func.nqName;
        let call: DG.FuncCall;
        try {
          call = func.prepare();
        }
        catch {
          unprepared++;
          continue;
        }
        const scope = Scope.root('func-form sweep');
        try {
          const form = scope.build(() => funcForm(call, {editors: 'auto'}));
          for (const p of form.problems.peek())
            problems.push(`${current}.${p.name}: ${p.message}`);
          let timer: ReturnType<typeof setTimeout> | undefined;
          const settled = await Promise.race([form.settled.then(() => true),
            new Promise<boolean>((resolve) => timer = setTimeout(() => resolve(false), SETTLE_MS))]);
          clearTimeout(timer);
          if (!settled)
            unsettled++;
        }
        catch (e) {
          thrown.push(`${current}: ${e instanceof Error ? e.message : e}`);
        }
        finally {
          scope.dispose();
        }
      }
    } finally {
      console.error = consoleError;
      window.removeEventListener('error', onError);
      window.removeEventListener('unhandledrejection', onRejection);
      grok.shell.closeTable(named);
      grok.shell.closeTable(nameless);
    }
    console.log(`funcForm sweep: ${funcs.length} functions, ${unprepared} could not be prepared, ` +
      `${unsettled} did not settle in ${SETTLE_MS} ms, ${problems.length} fields shown as errors` +
      (problems.length > 0 ? `: ${problems.join('; ')}` : ''));
    expect(thrown.length, 0, `funcForm threw for ${thrown.length}: ${thrown.slice(0, 10).join('; ')}`);
    expect(errors.length, 0, `console errors for ${errors.length}: ${errors.slice(0, 10).join('; ')}`);
  }, {timeout: 900000});
});
