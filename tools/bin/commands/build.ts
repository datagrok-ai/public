import fs from 'fs';
import path from 'path';
import {promisify} from 'util';
import {exec} from 'child_process';
import * as readline from 'readline';
import * as utils from '../utils/utils';
import * as color from '../utils/color-utils';

const execAsync = promisify(exec);

interface BuildArgs {
  _: string[];
  recursive?: boolean;
  silent?: boolean;
  filter?: string;
  verbose?: boolean;
  'no-incremental'?: boolean;
  parallel?: number;
}

interface BuildResult {
  name: string;
  version: string;
  buildTime: string;
  bundleSize: string;
  success: boolean;
}

interface PackageInfo {
  dir: string;
  name: string;
  friendlyName: string;
  version: string;
  packageJson: any;
}

export async function build(args: BuildArgs): Promise<boolean> {
  if (args.verbose)
    color.setVerbose(true);

  const buildCmd = args['no-incremental'] ? 'npm run build' : 'npm run build -- --env incremental';

  if (args.recursive)
    return await buildRecursive(process.cwd(), args, buildCmd);
  else
    return await buildSingle(process.cwd(), buildCmd);
}

async function buildSingle(dir: string, buildCmd: string): Promise<boolean> {
  if (!utils.isPackageDir(dir)) {
    color.error('Not a package directory (no package.json found)');
    return false;
  }

  const packageJson = JSON.parse(fs.readFileSync(path.join(dir, 'package.json'), 'utf-8'));
  const name = packageJson.friendlyName || packageJson.name;
  console.log(`Building ${name}...`);

  try {
    await utils.runScript('npm install', dir, color.isVerbose());
    await utils.runScript(buildCmd, dir, color.isVerbose());
    color.success(`Successfully built ${name}`);
    return true;
  }
  catch (error: any) {
    color.error(`Failed to build ${name}`);
    if (error.message)
      color.error(error.message);
    return false;
  }
}

async function buildRecursive(baseDir: string, args: BuildArgs, buildCmd: string): Promise<boolean> {
  const packages = discoverPackages(baseDir);
  if (packages.length === 0) {
    color.warn('No packages found in the current directory');
    return false;
  }

  const filtered = args.filter ? applyFilter(packages, args.filter) : packages;
  if (filtered.length === 0) {
    color.warn('No packages match the filter');
    return false;
  }

  console.log(`Found ${filtered.length} package(s): ${filtered.map((p) => p.friendlyName).join(', ')}`);

  if (!args.silent) {
    const confirmed = await confirm(`\nBuild ${filtered.length} package(s)?`);
    if (!confirmed) {
      console.log('Aborted.');
      return false;
    }
  }

  const maxParallel = args.parallel || 4;
  const results = await buildParallel(filtered, buildCmd, maxParallel);
  return results.every((r) => r.success);
}

function discoverPackages(baseDir: string): PackageInfo[] {
  const entries = fs.readdirSync(baseDir);
  const packages: PackageInfo[] = [];

  for (const entry of entries) {
    if (entry.startsWith('.') || entry === 'node_modules')
      continue;

    const dir = path.join(baseDir, entry);
    try {
      if (!fs.statSync(dir).isDirectory())
        continue;
    }
    catch (_) {
      continue;
    }

    const packageJsonPath = path.join(dir, 'package.json');
    if (!fs.existsSync(packageJsonPath))
      continue;

    try {
      const packageJson = JSON.parse(fs.readFileSync(packageJsonPath, 'utf-8'));
      packages.push({
        dir,
        name: packageJson.name || entry,
        friendlyName: packageJson.friendlyName || packageJson.name || entry,
        version: packageJson.version || '0.0.0',
        packageJson,
      });
    }
    catch (_) {
      continue;
    }
  }

  return packages.sort((a, b) => a.friendlyName.localeCompare(b.friendlyName));
}

function applyFilter(packages: PackageInfo[], filterStr: string): PackageInfo[] {
  const conditions = filterStr.split('&&').map((s) => s.trim());
  const parsedConditions = conditions.map((cond) => {
    const colonIdx = cond.indexOf(':');
    if (colonIdx === -1)
      return {field: cond, pattern: new RegExp('.')};

    const field = cond.substring(0, colonIdx).trim();
    const pattern = new RegExp(cond.substring(colonIdx + 1).trim());
    return {field, pattern};
  });

  return packages.filter((pkg) => {
    for (const cond of parsedConditions) {
      const value = getNestedValue(pkg.packageJson, cond.field);
      if (value === undefined || !cond.pattern.test(String(value)))
        return false;
    }
    return true;
  });
}

function getNestedValue(obj: any, path: string): any {
  const parts = path.split('.');
  let current = obj;
  for (const part of parts) {
    if (current == null || typeof current !== 'object')
      return undefined;
    current = current[part];
  }
  return current;
}

async function buildParallel(packages: PackageInfo[], buildCmd: string, maxParallel: number): Promise<BuildResult[]> {
  const results: BuildResult[] = [];

  const headers = ['Plugin', 'Version', 'Build time', 'Bundle size'];
  const widths = [
    Math.max(headers[0].length, ...packages.map((p) => p.friendlyName.length)),
    Math.max(headers[1].length, ...packages.map((p) => p.version.length)),
    Math.max(headers[2].length, 10),
    Math.max(headers[3].length, 40),
  ];
  const pad = (s: string, w: number) => s + ' '.repeat(Math.max(0, w - s.length));

  console.log(`\nBuilding with ${maxParallel} parallel job(s)...`);
  console.log(headers.map((h, i) => pad(h, widths[i])).join(' | '));
  console.log(widths.map((w) => '-'.repeat(w)).join('-+-'));

  const buildOne = async (pkg: PackageInfo): Promise<BuildResult> => {
    const start = Date.now();
    let success = true;
    let buildTime = '';
    let bundleSize = '';

    try {
      await execAsync('npm install', {cwd: pkg.dir, maxBuffer: 10 * 1024 * 1024});
      await execAsync(buildCmd, {cwd: pkg.dir, maxBuffer: 10 * 1024 * 1024});
      const elapsed = (Date.now() - start) / 1000;
      buildTime = `${elapsed.toFixed(1)}s`;
      bundleSize = getBundleSize(pkg.dir);
    }
    catch (error: any) {
      success = false;
      buildTime = 'Error';
      const raw = (error.stderr || error.stdout || error.message || 'Unknown error').trim();
      bundleSize = raw.replace(/\r?\n/g, ' ').replace(/\s+/g, ' ').substring(0, 40);
    }

    const result: BuildResult = {
      name: pkg.friendlyName,
      version: pkg.version,
      buildTime,
      bundleSize,
      success,
    };
    results.push(result);

    const cells = [result.name, result.version, result.buildTime, result.bundleSize];
    const line = cells.map((cell, j) => pad(cell, widths[j])).join(' | ');
    if (success)
      color.info(line);
    else
      color.error(line);

    return result;
  };

  let idx = 0;
  const next = async (): Promise<void> => {
    while (idx < packages.length) {
      const pkg = packages[idx++];
      await buildOne(pkg);
    }
  };
  const workers = Array.from({length: Math.min(maxParallel, packages.length)}, () => next());
  await Promise.all(workers);

  const succeeded = results.filter((r) => r.success).length;
  const failed = results.length - succeeded;
  console.log('');
  if (failed === 0)
    color.success(`All ${results.length} package(s) built successfully`);
  else
    color.warn(`${succeeded} succeeded, ${failed} failed`);

  return results;
}

function getBundleSize(dir: string): string {
  const bundlePath = path.join(dir, 'dist', 'package.js');
  try {
    const stats = fs.statSync(bundlePath);
    return `${(stats.size / 1024).toFixed(1)} KB`;
  }
  catch (_) {
    return 'N/A';
  }
}


function confirm(message: string): Promise<boolean> {
  const rl = readline.createInterface({input: process.stdin, output: process.stdout});
  return new Promise((resolve) => {
    rl.question(`${message} [Y/n] `, (answer) => {
      rl.close();
      resolve(answer.trim().toLowerCase() !== 'n');
    });
  });
}
