import fs from 'fs';
import os from 'os';
import path from 'path';
import {exec} from 'child_process';
import yaml from 'js-yaml';
import * as utils from '../utils/utils';
import * as color from '../utils/color-utils';
import {processPackage} from './publish';

const grokDir = path.join(os.homedir(), '.grok');
const confPath = path.join(grokDir, 'config.yaml');

interface RunArgs {
  _: string[];
  key?: string;
  release?: boolean;
  verbose?: boolean;
}

function getWebUrl(apiUrl: string): string {
  const u = new URL(apiUrl);
  u.pathname = u.pathname.replace(/\/api\/?$/, '') || '/';
  return u.toString().replace(/\/$/, '');
}

function openBrowser(url: string): void {
  let command: string;
  switch (process.platform) {
    case 'darwin':  command = `open "${url}"`; break;
    case 'win32':   command = `start "" "${url}"`; break;
    default:        command = `xdg-open "${url}"`;
  }
  exec(command, (err) => {
    if (err)
      color.warn(`Could not open browser: ${err.message}`);
  });
}

const MISSING_MODULE_PATTERNS = ['cannot find module', 'module not found', 'can\'t resolve'];

async function buildPackage(dir: string): Promise<boolean> {
  const buildCmd = 'npm run build -- --env incremental';
  const packageJson = JSON.parse(fs.readFileSync(path.join(dir, 'package.json'), 'utf-8'));
  const name = packageJson.friendlyName || packageJson.name;
  console.log(`Building ${name}...`);
  try {
    await utils.runScript(buildCmd, dir, color.isVerbose());
    color.success(`Successfully built ${name}`);
    return true;
  }
  catch (error: any) {
    const msg: string = (error?.message ?? '').toLowerCase();
    if (MISSING_MODULE_PATTERNS.some((p) => msg.includes(p))) {
      color.warn('Missing modules detected, running npm install...');
      try {
        await utils.runScript('npm install', dir, color.isVerbose());
        await utils.runScript(buildCmd, dir, color.isVerbose());
        color.success(`Successfully built ${name}`);
        return true;
      }
      catch (retryError: any) {
        color.error(`Failed to build ${name}`);
        if (retryError.message)
          color.error(retryError.message);
        return false;
      }
    }
    color.error(`Failed to build ${name}`);
    if (error.message)
      color.error(error.message);
    return false;
  }
}

export async function run(args: RunArgs): Promise<boolean> {
  color.setVerbose(args.verbose || false);

  // Step 1: Build (skip npm install to preserve npm link; retry with npm install only if needed)
  const built = await buildPackage(process.cwd());
  if (!built)
    return false;

  // Step 2: Resolve server URL and key
  if (!fs.existsSync(confPath)) {
    color.error(`Config not found at ${confPath}. Run \`grok config\` first.`);
    return false;
  }
  const config = yaml.load(fs.readFileSync(confPath, {encoding: 'utf-8'})) as utils.Config;
  const urls = utils.mapURL(config);

  let host = config.default;
  if (args['_'].length === 2) host = args['_'][1];

  let key = '';
  let url = '';
  let registry: string | undefined;

  try {
    url = new URL(host).href;
    if (url.endsWith('/')) url = url.slice(0, -1);
    if (url in urls) {
      const alias = urls[url];
      key = config['servers'][alias]['key'];
      registry = config['servers'][alias]['registry'];
    }
  }
  catch (error) {
    if (!(host in config.servers)) {
      color.error(`Unknown server alias. Please add it to ${confPath}`);
      return false;
    }
    url = config['servers'][host]['url'];
    key = config['servers'][host]['key'];
    registry = config['servers'][host]['registry'];
  }

  if (args.key) key = args.key;
  if (key === '') {
    color.warn('Please provide the key with `--key` option or add it by running `grok config`');
    return false;
  }

  const packDir = path.join(process.cwd(), 'package.json');
  if (!fs.existsSync(packDir)) {
    color.error('`package.json` doesn\'t exist');
    return false;
  }
  const packageName = JSON.parse(fs.readFileSync(packDir, {encoding: 'utf-8'})).name;

  // Step 3: Publish
  process.env.NODE_TLS_REJECT_UNAUTHORIZED = '0';
  const code = await processPackage(!args.release, false, url, key, packageName, false, undefined, host, registry);
  if (code !== 0)
    return false;

  // Step 4: Open browser
  const webUrl = getWebUrl(url);
  color.success(`Opening ${webUrl}`);
  openBrowser(webUrl);

  return true;
}
