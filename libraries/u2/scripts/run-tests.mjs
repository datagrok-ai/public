/* `node --test` over ./tests/*.test.js, listed here because cmd.exe does not expand a glob. The tests
   import `../src/**.js`, which `npm test` compiles into `.test-build/` rather than next to the sources;
   a resolve hook loads them from there. With `--dist-of <lib>` it first fails when <lib>/dist is missing
   or older than <lib>/src (ems runs against u2/dist). Other arguments go to `node --test`. A u2 usage error
   throws (`U2_STRICT=1`) unless `U2_STRICT` says otherwise (`report` lists them all). The DOM is happy-dom
   (tests/happy-dom.mjs); `U2_DOM=shim` runs the hand-written tests/dom-shim.js instead. */
import {existsSync, readdirSync, realpathSync, statSync} from 'node:fs';
import {join} from 'node:path';
import {pathToFileURL} from 'node:url';
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

if (!existsSync(join('.test-build', 'src'))) {
  console.error('.test-build/src is missing; the tests import that build. Run `npm test`, which compiles it first.');
  process.exit(1);
}

const root = pathToFileURL(realpathSync('.')).href + '/';
const hook = 'data:text/javascript,' + encodeURIComponent(`
const SRC = ${JSON.stringify(root + 'src/')};
const BUILT = ${JSON.stringify(root + '.test-build/src/')};
export async function resolve(specifier, context, next) {
  if ((specifier.startsWith('.') || specifier.startsWith('file:')) && context.parentURL?.startsWith('file:')) {
    const url = new URL(specifier, context.parentURL);
    if (url.href.startsWith(SRC) && url.pathname.endsWith('.js'))
      return next(BUILT + url.href.slice(SRC.length), context);
  }
  return next(specifier, context);
}
`);
const register = 'data:text/javascript,' +
  encodeURIComponent(`import {register} from 'node:module';\nregister(${JSON.stringify(hook)});\n`);

const dom = process.env.U2_DOM === 'shim' ? [] : ['--import', new URL('../tests/happy-dom.mjs', import.meta.url).href];
const files = readdirSync('tests').filter((f) => f.endsWith('.test.js')).sort().map((f) => `tests/${f}`);
const r = spawnSync(process.execPath, [...dom, '--import', register, '--test', ...args, ...files],
  {stdio: 'inherit', env: {...process.env, U2_STRICT: process.env.U2_STRICT ?? '1'}});
process.exit(r.status ?? 1);
