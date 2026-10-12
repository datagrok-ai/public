/** Progress of long batches on the task bar: "what: 45,000 of 100,000". */
import * as DG from 'datagrok-api/dg';
import type {Progress} from './pool';

/** Report `done of total` on `pi`. */
export function taskProgress(pi: DG.TaskBarProgressIndicator, what: string): Progress {
  return (done, total) => pi.update(Math.round(100 * done / total),
    `${what}: ${done.toLocaleString()} of ${total.toLocaleString()}`);
}

/** Run `body` with its own task bar progress when `total` is large (functions that
 * other dialogs call, e.g. engines of Bio's numbering and MSA dialogs). */
export async function withProgress<T>(what: string, total: number, body: (progress?: Progress) => Promise<T>):
  Promise<T> {
  if (total < 2000) return body();
  const pi = DG.TaskBarProgressIndicator.create(`${what}...`);
  try {
    return await body(taskProgress(pi, what));
  } finally {
    pi.close();
  }
}
