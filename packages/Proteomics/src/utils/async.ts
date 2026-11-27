/** Polls `check` every `interval` ms until it returns true, rejecting with
 * `error` after `wait` ms. Production counterpart of the test library's
 * awaitCheck (which must not be imported outside tests). */
export function awaitCheck(
  check: () => boolean, error: string = 'Timeout exceeded', wait: number = 500, interval: number = 50,
): Promise<void> {
  return new Promise((resolve, reject) => {
    const started = Date.now();
    const tick = (): void => {
      let ok = false;
      try {ok = check();} catch {ok = false;}
      if (ok) {resolve(); return;}
      if (Date.now() - started >= wait) {reject(new Error(error)); return;}
      setTimeout(tick, interval);
    };
    tick();
  });
}
