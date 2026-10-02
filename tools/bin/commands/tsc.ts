import {spawnSync} from 'child_process';
import {toolchain, copyLibraryAssets, ensureLibraryExports} from '../utils/toolchain';

/**
 * `grok tsc ...` runs the workspace TypeScript compiler (the one `@datagrok/build-config` pins) with the
 * given arguments, so package scripts never depend on a locally installed `tsc`. After an emitting run
 * in a library, the assets its sources import are mirrored into dist/.
 */
export async function tsc(args: {_: string[]}): Promise<boolean> {
  const cwd = process.cwd();
  const extra = process.argv.slice(process.argv.indexOf('tsc') + 1);
  const emitting = !extra.includes('--noEmit');
  if (emitting)
    ensureLibraryExports(cwd);
  const r = spawnSync(process.execPath, [toolchain(cwd).tsc, ...extra], {stdio: 'inherit', cwd});
  if (r.status !== 0)
    return false;
  if (emitting)
    copyLibraryAssets(cwd);
  return true;
}
