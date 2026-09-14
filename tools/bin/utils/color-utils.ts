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

/**
 * One step of a multi-step command: prints its label, ticks a spinner while [action] runs, and
 * replaces the line with the outcome. A spinner needs a terminal to erase lines, so a CI log
 * (or a redirected stdout) gets one plain line per step instead.
 */
export async function step<T>(label: string, action: () => Promise<T>): Promise<T> {
  const tty = process.stdout.isTTY === true;
  const frames = ['-', '\\', '|', '/'];
  let frame = 0;
  const draw = () => process.stdout.write(`\r  ${frames[frame++ % frames.length]} ${label}   `);
  if (!tty)
    console.log(`  ${label}...`);
  const timer = tty ? setInterval(draw, 120) : null;
  if (tty)
    draw();
  const finish = (mark: string, color: string, text: string) => {
    if (timer)
      clearInterval(timer);
    if (tty)
      process.stdout.write(`\r\x1b[2K  \x1b[${color}m${mark}\x1b[0m ${text}\n`);
    else if (mark !== '+')
      console.log(`  ${mark} ${text}`);
  };
  try {
    const result = await action();
    finish('+', '32', label);
    return result;
  } catch (e) {
    finish('x', '31', label);
    throw e;
  }
}
