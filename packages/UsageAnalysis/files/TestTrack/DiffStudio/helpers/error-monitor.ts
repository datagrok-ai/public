import { Page } from '@playwright/test';

export interface CapturedError { source: string; text: string; }

const priorMonitors = new WeakMap<Page, Array<{ event: string; fn: (...args: any[]) => void }>>();

export function attachErrorMonitor(page: Page): {
  errors: CapturedError[];
  assertNone: () => void;
} {
  const errors: CapturedError[] = [];

  const prev = priorMonitors.get(page);
  if (prev) {
    for (const { event, fn } of prev) page.off(event as any, fn as any);
  }
  const registered: Array<{ event: string; fn: (...args: any[]) => void }> = [];

  const isNoise = (text: string): boolean => {

    if (text.includes('Failed to load resource') && text.includes('404')) return true;
    if (text.includes('Failed to load resource') && text.includes('502')) return true;
    if (text.includes('Failed to load resource') && text.includes('sourcemap')) return true;

    if (/^Failed to load resource: the server responded with a status of \d+ \(.*\)$/.test(text)) {
      return true;
    }

    if (text.includes('Jiraconnect')) return true;
    if (text.includes('Error fetching projects')) return true;

    if (text.includes('Failed to fetch')) return true;

    if (text.includes('WebSocket connection') && text.includes('failed')) return true;

    if (text.includes('packages/$sdk/lib/')) return true;
    if (text.includes('_FutureListener.handleError')) return true;
    if (text.includes('_propagateToListeners')) return true;

    if (text.includes('Plates:GetPlates')) return true;
    if (text.includes("Failed host lookup: 'database'")) return true;
    if (text.includes('SocketException: Failed host lookup')) return true;

    if (text.includes('connectors/connections/System.AppData/file/DiffStudio')) return true;

    if (text.includes('WebAssembly compilation aborted')) return true;
    if (text.includes('Response body loading was aborted')) return true;

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

export async function snapshotBalloonErrors(page: Page): Promise<string[]> {
  return await page.evaluate(() =>
    Array.from(document.querySelectorAll('.d4-balloon-error'))
      .map(e => (e.textContent ?? '').trim())
      .filter(t => t.length > 0));
}
