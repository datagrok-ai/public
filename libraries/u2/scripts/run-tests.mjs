/* `node --test` over ./tests/*.test.js, listed here because cmd.exe does not expand a glob. With
   `--dist-of <lib>` it first fails when <lib>/dist is missing or older than <lib>/src (ems runs
   against u2/dist). Other arguments go to `node --test`. */
import {existsSync, readdirSync, statSync} from 'node:fs';
import {join} from 'node:path';
import {spawnSync} from 'node:child_process';

const args = process.argv.slice(2);
const at = args.indexOf('--dist-of');
const lib = at >= 0 ? args.splice(at, 2)[1] : null;

function newestSource(dir) {
  let t = 0;
  for (const e of readdirSync(dir, {withFileTypes: true})) {
    const p = join(dir, e.name);
    if (e.isDirectory())
      t = Math.max(t, newestSource(p));
    else if (e.name.endsWith('.ts') && !e.name.endsWith('.d.ts'))
      t = Math.max(t, statSync(p).mtimeMs);
  }
  return t;
}

if (lib) {
  const built = join(lib, 'dist', 'src', 'index.js');
  if (!existsSync(built) || statSync(built).mtimeMs < newestSource(join(lib, 'src'))) {
    console.error(`${built} is ${existsSync(built) ? 'older than the sources in ' + join(lib, 'src') : 'missing'}; ` +
      `these tests import that build. Run \`pnpm -C ${lib} build\` (or \`pnpm turbo run test\`) first.`);
    process.exit(1);
  }
}

const files = readdirSync('tests').filter((f) => f.endsWith('.test.js')).sort().map((f) => `tests/${f}`);
const r = spawnSync(process.execPath, ['--test', ...args, ...files], {stdio: 'inherit'});
process.exit(r.status ?? 1);
