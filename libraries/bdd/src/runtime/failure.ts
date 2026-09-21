/* What a failed step reports: the feature line, the step as written, and the reason in one plain
   sentence — never a matcher diff, an in-page stack or a harness line. Pure, so it is unit-tested. */

const PROGRAMMING_ERRORS = ['TypeError', 'ReferenceError', 'RangeError', 'SyntaxError'];
const INTERNAL_FRAME = /node_modules[\\/]+(?:playwright|@playwright)|UtilityScript|eval at evaluate|node:internal|node:async_hooks/;
const ANSI = new RegExp(String.fromCharCode(27) + '\\[[0-9;]*m', 'g');

/** The failure of one Gherkin step. `message` is the report; the stack carries frames only for a
 * programming error in a binding, where the frame is the information. */
export class StepFailure extends Error {
  /** `frame` is the feature line as a stack frame (`<abs path>:<line>:1`), so a reporter that
   * prints a snippet prints the Gherkin. */
  constructor(readonly at: string, readonly step: string, readonly reason: string, cause?: unknown, readonly frame = '') {
    super(`${at}\n  ${step}\n\n${reason}`);
    this.name = 'StepFailure';
    this.stack = this.message + (frames(cause) || (frame ? `\n    at ${frame}` : ''));
  }
}

/** The reason in an error, as a sentence: Playwright's API prefix ("locator.evaluate: Error: "),
 * the in-page stack of an evaluate, the millisecond timeout phrasing and the selector it was
 * waiting for (the step names the phrase) go; a matcher's own explanation (expected/received,
 * what the element resolved to, why the action did not go through) stays. */
export function reasonOf(e: unknown): string {
  const text = (e instanceof Error ? e.message : String(e)).replace(ANSI, '');
  const lines = text.split('\n').filter((l) => !/^\s+at\s/.test(l) && !/^\s+- waiting for (?:locator|getBy|internal:)/.test(l));
  while (lines.length > 0 && lines[lines.length - 1].trim() === '')
    lines.pop();
  if (lines[lines.length - 1]?.trim() === 'Call log:')
    lines.pop();
  let s = lines.join('\n').trim();
  s = s.replace(/^(?:page|locator|frame|mouse|keyboard|expect)\.\w+: (?:Error: )?/, '');
  s = s.replace(/^Timeout (\d+)ms exceeded\.?/, (_, ms) => `timed out after ${Number(ms) / 1000} s`);
  return s;
}

/** Whether the error is Playwright giving up on an element — the case where what the page shows
 * instead is the missing half of the report. */
export function isWaitFailure(e: unknown): boolean {
  // Playwright colours its matcher output, and a check that names itself puts its message before
  // the "expect(locator)… failed" line
  const text = (e instanceof Error ? e.message : String(e)).replace(/\x1b\[[0-9;]*m/g, '');
  return (e instanceof Error && e.name === 'TimeoutError') || /Timeout \d+ms exceeded|(^|\n)expect\(/.test(text);
}

export function failure(at: string, step: string, e: unknown, shown = '', frame = ''): StepFailure {
  if (e instanceof StepFailure)
    return e;
  const reason = reasonOf(e);
  return new StepFailure(at, step, shown ? `${reason}\n${shown}` : reason, e, frame);
}

/** The journey's verdict: every failed scenario with its step report; the first failure's feature
 * line is the only frame. */
export function journeyFailure(failed: {name: string; error: unknown}[], scenarios: number): Error {
  const list = failed.map((f) => `${f.name}\n${indent(f.error instanceof StepFailure ? f.error.message : reasonOf(f.error))}`);
  const e = new Error(`${failed.length} of ${scenarios} scenarios failed\n\n${list.join('\n\n')}`);
  const first = failed.map((f) => f.error).find((x): x is StepFailure => x instanceof StepFailure && x.frame !== '');
  e.stack = e.message + (first ? `\n    at ${first.frame}` : '');
  return e;
}

function indent(s: string): string {
  return s.split('\n').map((l) => l.length > 0 ? `  ${l}` : l).join('\n');
}

function frames(e: unknown): string {
  if (!(e instanceof Error) || !PROGRAMMING_ERRORS.includes(e.name) || !e.stack)
    return '';
  const kept = e.stack.split('\n').filter((l) => /^\s+at\s/.test(l) && !INTERNAL_FRAME.test(l));
  return kept.length > 0 ? `\n${kept.join('\n')}` : '';
}
