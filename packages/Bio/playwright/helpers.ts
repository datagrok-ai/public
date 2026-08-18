import * as fs from 'fs';
import * as os from 'os';
import * as path from 'path';

// Monomer-library selection lives in per-user server settings, so the specs that toggle libraries
// cannot run next to each other: one spec's excluded library breaks the other's detection or
// checkbox state. Playwright serialises only within a file, so take a lock across worker processes.
const LOCK_DIR = path.join(os.tmpdir(), 'dg-bio-monomer-lib.lock');
const STALE_MS = 5 * 60_000;
const WAIT_MS = 4 * 60_000;

export async function acquireMonomerLibLock(): Promise<void> {
  const deadline = Date.now() + WAIT_MS;
  for (;;) {
    try {
      fs.mkdirSync(LOCK_DIR);
      return;
    }
    catch (e) {
      let ageMs = 0;
      try { ageMs = Date.now() - fs.statSync(LOCK_DIR).mtimeMs; }
      catch (_) { continue; }
      // a lock left behind by a killed run must not wedge every later run
      if (ageMs > STALE_MS) {
        fs.rmSync(LOCK_DIR, {recursive: true, force: true});
        continue;
      }
      if (Date.now() > deadline)
        throw new Error(`monomer-library lock held for over ${WAIT_MS / 1000}s by another spec`);
      await new Promise((r) => setTimeout(r, 1000));
    }
  }
}

export function releaseMonomerLibLock(): void {
  fs.rmSync(LOCK_DIR, {recursive: true, force: true});
}
