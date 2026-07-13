import { test } from '@playwright/test';

export interface StepError { step: string; error: string; }

/** Returns a per-test soft-step collector that captures failures and asserts them at the end. */
export function createSoftStepCollector() {
  const errors: StepError[] = [];

  const softStep = async (name: string, fn: () => Promise<void>): Promise<void> => {
    try {
      await test.step(name, fn);
    } catch (e: any) {
      errors.push({ step: name, error: e?.message ?? String(e) });
      console.error(`[STEP FAILED] ${name}: ${e?.message ?? e}`);
    }
  };

  const assertAllPassed = (): void => {
    if (errors.length > 0) {
      const summary = errors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
      throw new Error(`${errors.length} step(s) failed:\n${summary}`);
    }
  };

  return { softStep, assertAllPassed };
}
