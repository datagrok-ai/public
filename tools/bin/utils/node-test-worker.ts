/* eslint-disable max-len */
/**
 * Standalone worker for the `grok test` Node (browserless) pass.
 *
 * Spawned by node-test-runner.ts with cwd = the package directory:
 * 1. loads the js-api Node runtime (`datagrok-api/datagrok`) from the package's own node_modules,
 * 2. resolves the published package version the platform would serve and evaluates its
 *    `dist/package-test.js` bundle (webpack var-library `<name>_test`),
 * 3. runs init-tagged functions of the test bundle, then calls its exported `testNode()`
 *    which executes only tests marked `{node: true}`,
 * 4. counts the tests the browser pass still has to run, and writes a JSON report to --out.
 *
 * Exit codes: 0 = pass executed (test failures are conveyed in the report),
 * 3 = Node pass unsupported for this package (caller falls back to the browser for everything),
 * anything else = infrastructure error.
 */
import * as vm from 'vm';
import * as fs from 'fs';
import * as path from 'path';

const UNSUPPORTED_EXIT_CODE = 3;
const RESULT_COLUMNS = ['date', 'category', 'name', 'success', 'result', 'ms', 'skipped', 'logs', 'owner', 'package', 'widgetsDifference', 'flaking'];

function getArg(name: string): string | undefined {
  const prefix = `--${name}=`;
  const arg = process.argv.find((a) => a.startsWith(prefix));
  return arg ? arg.slice(prefix.length) : undefined;
}

function hasFlag(name: string): boolean {
  return process.argv.includes(`--${name}`);
}

function unsupported(reason: string): never {
  console.log(`Node pass unavailable: ${reason}`);
  process.exit(UNSUPPORTED_EXIT_CODE);
}

// A crashing async callback (e.g. a socket handler inside the Dart runtime) must not
// kill the whole pass — log it and let the per-test timeout fail the affected test.
process.on('uncaughtException', (err) => {
  console.error('Uncaught exception:', err?.stack ?? err);
});
process.on('unhandledRejection', (reason) => {
  console.error('Unhandled rejection:', reason);
});

async function main(): Promise<void> {
  const apiUrl = process.env.DG_API_URL!;
  const token = process.env.DG_API_TOKEN!;
  const packageName = process.env.TARGET_PACKAGE!;
  const outPath = getArg('out')!;
  const categoryArg = getArg('category');
  const testArg = getArg('test');
  const stressTest = hasFlag('stress-test');
  const benchmark = hasFlag('benchmark');
  const verbose = hasFlag('verbose');

  let dgEntry: string;
  try {
    dgEntry = require.resolve('datagrok-api/datagrok', {paths: [process.cwd()]});
  } catch {
    unsupported('the installed datagrok-api has no Node runtime (datagrok-api/datagrok)');
  }
  const dgApi = require(dgEntry!);
  if (typeof dgApi.startDatagrok !== 'function')
    unsupported('the installed datagrok-api has no startDatagrok()');

  await dgApi.startDatagrok({apiUrl, apiToken: token});
  const g: any = globalThis;
  const grok = g.grok;
  const DG = g.DG;

  // Package bundles map webpack externals to browser globals (rxjs, rxjs.operators,
  // dayjs, utc, wu, $). The Node runtime guarantees DG/grok/ui; provide the rest
  // from datagrok-api's own dependencies so bundle evaluation doesn't throw.
  const req = (id: string) => require(require.resolve(id, {paths: [path.dirname(dgEntry!), process.cwd()]}));
  const ensureGlobal = (name: string, resolve: () => any) => {
    if (g[name] == null)
      try {
        g[name] = resolve();
      } catch {}
  };
  ensureGlobal('rxjs', () => req('rxjs'));
  if (g.rxjs != null && g.rxjs.operators == null)
    try {
      g.rxjs = Object.assign({}, g.rxjs, {operators: req('rxjs/operators')});
    } catch {}
  ensureGlobal('dayjs', () => req('dayjs'));
  ensureGlobal('utc', () => req('dayjs/plugin/utc'));
  ensureGlobal('wu', () => req('wu'));
  ensureGlobal('$', () => req('cash-dom'));

  // Resolve the published version the platform would serve: the caller's debug version
  // if one exists, else the current one (mirrors js-api's Node loadPackage()).
  const headers = {Authorization: token};
  const listResp = await fetch(`${apiUrl}/packages?text=${encodeURIComponent(packageName)}`, {headers});
  if (!listResp.ok)
    throw new Error(`Failed to list packages: HTTP ${listResp.status}`);
  const packages: any[] = await listResp.json();
  const pkgInfo = packages.find((p) => p.name?.toLowerCase() === packageName.toLowerCase());
  if (!pkgInfo)
    unsupported(`package "${packageName}" is not published on ${apiUrl}`);
  const versions: any[] = pkgInfo.publishedVersions ?? [];
  const login = (await grok.dapi.users.current()).login;
  const cur =
    versions.find((v) => v.debug && v.version === login) ??
    versions.find((v) => v.isCurrent) ??
    versions.find((v) => v.isLatest) ??
    versions[0];
  if (!cur)
    unsupported(`package "${packageName}" has no published versions`);

  const webRoot = `${apiUrl}/packages/published/files/${pkgInfo.name}/${cur.version}/${cur.buildHash}/${cur.buildNumber}/`;
  const jsUrl = `${webRoot}dist/package-test.js`;
  const codeResp = await fetch(jsUrl, {headers});
  if (!codeResp.ok)
    unsupported(`no test bundle at ${jsUrl} (HTTP ${codeResp.status})`);
  vm.runInThisContext(await codeResp.text(), {filename: jsUrl});

  const moduleName = `${packageName.toLowerCase()}_test`;
  const testModule = g[moduleName] ?? g[Object.keys(g).find((k) => k.endsWith('_test') && g[k]?.tests) ?? ''];
  if (!testModule?.tests)
    unsupported(`test bundle did not define a "${moduleName}" module`);
  if (typeof testModule.testNode !== 'function')
    unsupported('package-test.ts does not export testNode()');

  if (testModule._package != null) {
    testModule._package.webRoot = webRoot;
    testModule._package.name = pkgInfo.name;
    testModule._package.version = cur.version;
  }

  // The package must carry the resolved published-version id: functions are attached to
  // package versions, and runTests' initAutoTests filters `package.id = <id>` to register
  // the annotation-driven auto tests. The family-level entity id matches nothing.
  const pkg = await grok.dapi.packages.find(cur.id) ??
    await grok.dapi.packages.filter(`shortName = "${pkgInfo.name}"`).first();
  if (!pkg)
    unsupported(`package entity "${pkgInfo.name}" not found`);
  try {
    if (pkg.name !== pkgInfo.name)
      pkg.name = pkgInfo.name;
  } catch {}

  if (benchmark)
    DG.Test.isInBenchmark = true;

  const results: any[] = await testModule.testNode(pkg, {category: categoryArg, test: testArg, stressTest, verbose});

  // Count what the browser pass still has to run. testNode() ran initAutoTests(), so the
  // registry now also holds the annotation-driven auto/demo/detector tests.
  let browserTestsRemaining = 0;
  for (const [catName, cat] of Object.entries<any>(testModule.tests)) {
    if (categoryArg != null && !catName.toLowerCase().startsWith(categoryArg.toLowerCase().trim()))
      continue;
    for (const t of cat.tests ?? []) {
      if (testArg != null && t.name.toLowerCase() !== testArg.toLowerCase())
        continue;
      if (benchmark && !t.options?.benchmark)
        continue;
      if (stressTest && !t.options?.stressTest)
        continue;
      if (!(t.options?.node ?? cat.node ?? false))
        browserTestsRemaining++;
    }
  }

  // Shape the report like the browser pass does — same CSV columns in the same order,
  // so the two passes merge row-wise.
  const Papa = require('papaparse');
  const rows = results.map((r) => {
    const row: any = {};
    for (const c of RESULT_COLUMNS)
      row[c] = r[c] ?? (c === 'widgetsDifference' ? 0 : '');
    return row;
  });
  let failed = false;
  let passedAmount = 0; let skippedAmount = 0; let failedAmount = 0;
  let verbosePassed = ''; let verboseSkipped = ''; let verboseFailed = '';
  for (const r of results) {
    const line = `${r.category}: ${r.name} (${r.ms} ms) :  ${r.result}\n`;
    if (r.skipped) {
      skippedAmount++;
      verboseSkipped += line;
    } else if (r.success) {
      passedAmount++;
      verbosePassed += line;
    } else {
      failedAmount++;
      verboseFailed += line;
      failed = true;
    }
  }
  const report = {
    failed, verbosePassed, verboseSkipped, verboseFailed,
    passedAmount, skippedAmount, failedAmount,
    csv: rows.length ? Papa.unparse(rows, {columns: RESULT_COLUMNS}) : '',
    browserTestsRemaining,
    nodeTestsRun: results.length,
  };
  fs.writeFileSync(outPath, JSON.stringify(report), 'utf8');
  process.exit(0);
}

main().catch((e: any) => {
  console.error(e?.stack ?? e?.message ?? String(e));
  process.exit(2);
});
