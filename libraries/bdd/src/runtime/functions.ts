/* A package function called through the platform and its last result, kept in the page for
   the checks that follow (`bindings/platform/functions.ts`, and a package's own readings of it). */
import type {Page} from '@playwright/test';

declare const grok: any;

/** The call's arguments from a `| name | value |` table: `column:X` is column X of the current
 * table, `table` the current table, a number is a number, `true`/`false` a boolean, `""` empty. */
export function callFunction(page: Page, name: string, rows: string[][]): Promise<void> {
  return page.evaluate(async ([n, args]) => {
    const w = window as any;
    const params: Record<string, unknown> = {};
    for (const [key, raw] of args) {
      const df = grok.shell.t;
      if (raw === 'table')
        params[key] = df;
      else if (raw.startsWith('column:')) {
        const col = df?.col(raw.slice(7));
        if (!col)
          throw new Error(`no "${raw.slice(7)}" column in the current table; it has: ${df?.columns.names().join(', ') ?? 'no table'}`);
        params[key] = col;
      }
      else if (raw === 'true' || raw === 'false')
        params[key] = raw === 'true';
      else if (raw !== '' && !isNaN(Number(raw)))
        params[key] = Number(raw);
      else
        params[key] = raw;
    }
    let result: unknown;
    try {
      result = await grok.functions.call(n, params);
    }
    catch (e: any) {
      throw new Error(`${n} failed: ${String(e?.message ?? e).split(/\r?\n/)[0]}`);
    }
    w.__bddLastResult = {name: n, value: result};
  }, [name, rows] as [string, string[][]]);
}

/** A reading of the last result made in the page (the value is a live platform object):
 * `expression` sees `value` and `arg`, and must return something serializable. */
export function readResult(page: Page, expression: string, arg: unknown = null): Promise<any> {
  return page.evaluate(([e, a]) => {
    const w = window as any;
    if (!w.__bddLastResult)
      throw new Error('no function has been called in this scenario');
    return new Function('value', 'arg', `return (${e});`)(w.__bddLastResult.value, a);
  }, [expression, arg] as [string, unknown]);
}
