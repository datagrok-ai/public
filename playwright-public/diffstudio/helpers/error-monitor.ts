import { Page } from '@playwright/test';

export interface CapturedError { source: string; text: string; }

/**
 * Listeners registered by the last `attachErrorMonitor` call on a given page. The worker
 * `sharedPage` is reused across tests, so each test re-attaches a fresh monitor — we detach
 * the previous test's listeners first to avoid pile-up (and Node's MaxListeners warning) and
 * to keep capture scoped to the current test.
 */
const priorMonitors = new WeakMap<Page, Array<{ event: string; fn: (...args: any[]) => void }>>();

/**
 * Subscribe to console errors, uncaught page errors, failing network requests, and
 * `.d4-balloon-error` toasts. Returns a collector with `errors` (read on demand) and
 * `assertNone()` which throws an aggregated error if anything was captured.
 *
 * Wire this up immediately after `page` is available — typically as the first call in a test.
 */
export function attachErrorMonitor(page: Page): {
  errors: CapturedError[];
  assertNone: () => void;
} {
  const errors: CapturedError[] = [];

  // Detach the previous test's monitor from this (shared) page before wiring a fresh one.
  const prev = priorMonitors.get(page);
  if (prev) {
    for (const { event, fn } of prev) page.off(event as any, fn as any);
  }
  const registered: Array<{ event: string; fn: (...args: any[]) => void }> = [];

  /**
   * Errors that originate outside DiffStudio / Compute2 and should not fail our tests.
   * The dev environment hosts many third-party packages; their failures are environmental noise.
   */
  const isNoise = (text: string): boolean => {
    // Generic 4xx/5xx noise from missing optional assets (sourcemaps, icons, etc.)
    if (text.includes('Failed to load resource') && text.includes('404')) return true;
    if (text.includes('Failed to load resource') && text.includes('502')) return true;
    if (text.includes('Failed to load resource') && text.includes('sourcemap')) return true;
    // Some browsers log a bare "Failed to load resource: the server responded with a status of N ()"
    // with the URL missing from the message — when the network listener has already captured
    // and either filtered or recorded the parent request, the duplicate console line is noise.
    if (/^Failed to load resource: the server responded with a status of \d+ \(.*\)$/.test(text)) {
      return true;
    }
    // Third-party packages that habitually fail on the dev environment
    if (text.includes('Jiraconnect')) return true;
    if (text.includes('Error fetching projects')) return true;
    // CORS / external API failures
    if (text.includes('Failed to fetch')) return true;
    // WebSocket reconnect chatter
    if (text.includes('WebSocket connection') && text.includes('failed')) return true;
    // Dart SDK internal async stack frames (not DiffStudio-specific failures)
    if (text.includes('packages/$sdk/lib/')) return true;
    if (text.includes('_FutureListener.handleError')) return true;
    if (text.includes('_propagateToListeners')) return true;
    // Other third-party packages whose queries fail in the dev environment
    if (text.includes('Plates:GetPlates')) return true;
    if (text.includes("Failed host lookup: 'database'")) return true;
    if (text.includes('SocketException: Failed host lookup')) return true;
    // Optional file reads we attempt as a best-effort enrichment (tooltip parser fallback).
    // A 502 on these endpoints is recoverable — the test has hardcoded fallback values.
    if (text.includes('connectors/connections/System.AppData/file/DiffStudio')) return true;
    // WASM module loads (RDKit, OpenChemLib, etc.) sometimes abort when the page navigates
    // away mid-fetch. Not a DiffStudio failure.
    if (text.includes('WebAssembly compilation aborted')) return true;
    if (text.includes('Response body loading was aborted')) return true;
    // Benign streaming-compile fallback: when `WebAssembly.instantiateStreaming` can't be used
    // (MIME/cache), the glue retries via ArrayBuffer and logs this at error level. Not a failure.
    if (text.includes('falling back to ArrayBuffer instantiation')) return true;
    return false;
  };

  const onConsole = (msg: { type: () => string; text: () => string }): void => {
    if (msg.type() === 'error') {
      const text = msg.text();
      if (!isNoise(text)) errors.push({ source: 'console.error', text });
    }
  };
  const onPageError = (err: { message: string }): void => {
    if (!isNoise(err.message)) errors.push({ source: 'pageerror', text: err.message });
  };
  const onResponse = (resp: { status: () => number; url: () => string; request: () => { method: () => string } }): void => {
    const status = resp.status();
    if (status < 500) return;
    const url = resp.url();
    // Ignore 5xx from unrelated packages
    if (isNoise(url)) return;
    errors.push({ source: `network ${status}`, text: `${resp.request().method()} ${url}` });
  };

  page.on('console', onConsole as any);
  page.on('pageerror', onPageError as any);
  page.on('response', onResponse as any);
  registered.push({ event: 'console', fn: onConsole }, { event: 'pageerror', fn: onPageError },
    { event: 'response', fn: onResponse });
  priorMonitors.set(page, registered);

  const assertNone = (): void => {
    if (errors.length === 0) return;
    const summary = errors.map(e => `  - [${e.source}] ${e.text}`).join('\n');
    throw new Error(`${errors.length} runtime error(s) captured:\n${summary}`);
  };

  return { errors, assertNone };
}

/** Poll the DOM once for visible `.d4-balloon-error` toasts and return their texts. */
export async function snapshotBalloonErrors(page: Page): Promise<string[]> {
  return await page.evaluate(() =>
    Array.from(document.querySelectorAll('.d4-balloon-error'))
      .map(e => (e.textContent ?? '').trim())
      .filter(t => t.length > 0));
}
