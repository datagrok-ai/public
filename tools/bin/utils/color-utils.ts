export const error = (s: string) => console.log('\x1b[31m%s\x1b[0m', s);
export const info = (s: string) => console.log('\x1b[32m%s\x1b[0m', s);
export const warn = (s: string) => console.log('\x1b[33m%s\x1b[0m', s);

export const success = info;
export const fail = error;

let verbose = false;

export const setVerbose = (value: boolean) => verbose = value;
export const isVerbose = () => verbose;
export type LogType = 'error' | 'fail' | 'warn' | 'success' | 'info' | 'plain';

/** Logs a message only when verbose mode is enabled
 * @param {string} s - The message to log
 * @param {LogType} type - The type of the message, which determines its color. Defaults to 'plain'.
 */
export function log(s: string, type: LogType = 'plain'): void {
  if (!verbose)
    return;

  switch (type) {
  case 'fail':
  case 'error':
    error(s);
    break;
  case 'warn':
    warn(s);
    break;
  case 'success':
  case 'info':
    info(s);
    break;
  case 'plain':
  default:
    console.log(s);
    break;
  }
}
