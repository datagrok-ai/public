/* A package function called through the platform (`bindings/platform/functions.ts`). */
import type {Page} from '@playwright/test';

declare const grok: any;

export function callFunction(page: Page, name: string): Promise<void> {
  return page.evaluate(async (n) => {
    try {
      await grok.functions.call(n, {});
    }
    catch (e: any) {
      throw new Error(`${n} failed: ${String(e?.message ?? e).split(/\r?\n/)[0]}`);
    }
  }, name);
}
