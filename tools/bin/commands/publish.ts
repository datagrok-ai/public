// @ts-ignore
import archiver from 'archiver-promise';
import fs from 'fs';
// @ts-ignore
import fetch from 'node-fetch';
import os from 'os';
import path from 'path';
import walk from 'ignore-walk';
import yaml from 'js-yaml';
import {getDevKey, getToken} from '../utils/test-utils';
import * as utils from '../utils/utils';
import {Indexable} from '../utils/utils';
import {loadPackages} from '../utils/test-utils';
import * as color from '../utils/color-utils';
import {check} from './check';

const {exec} = require('child_process');

const grokDir = path.join(os.homedir(), '.grok');
const confPath = path.join(grokDir, 'config.yaml');
const confTemplateDir = path.join(path.dirname(path.dirname(__dirname)), 'config-template.yaml');
const confTemplate = yaml.load(fs.readFileSync(confTemplateDir, {encoding: 'utf-8'}));
let curDir = process.cwd();
const packDir = path.join(curDir, 'package.json');

const packageFiles = [
  'src/package.ts', 'src/detectors.ts', 'src/package.js', 'src/detectors.js',
  'src/package-test.ts', 'src/package-test.js', 'package.js', 'detectors.js',
];
let config: utils.Config;

export async function processPackage(debug: boolean, rebuild: boolean, host: string, devKey: string, packageName: any, dropDb: boolean, suffix?: string, hostAlias?: string) {
  // Get the server timestamps
  let timestamps: Indexable = {};
  let url = `${host}/packages/dev/${devKey}/${packageName}`;
  if (debug) {
    try {
      timestamps = await (await fetch(url + '/timestamps')).json();
      if (timestamps['#type'] === 'ApiError') {
        color.error(timestamps.message);
        return 1;
      }
    } catch (error) {
      if (utils.isConnectivityError(error))
        color.error(`Server is possibly offline: ${host}`);
      if (color.isVerbose())
        console.error(error);
      return 1;
    }
  }

  const zip = archiver('zip', {store: false});
  const chunks = [];
  zip.on('data', (chunk: any) => chunks.push(chunk));

  // Gather the files
  const localTimestamps: Indexable = {};
  const files = await walk({
    path: '.',
    ignoreFiles: ['.npmignore', '.gitignore', '.grokignore'],
    includeEmpty: false,
    follow: true,
  });

  const isWebpack = fs.existsSync('webpack.config.js');
  if (!rebuild && isWebpack) {
    if (fs.existsSync('dist/package.js')) {
      const distFiles = await walk({
        path: './dist',
        ignoreFiles: [],
        includeEmpty: false,
        follow: true,
      });
      distFiles.forEach((df: any) => {
        files.push(`dist/${df}`);
      });
    } else {
      color.warn('File `dist/package.js` not found. Building the package on the server side...\n' +
        'Next time, please build your package locally with Webpack beforehand\n' +
        'or run `grok publish` with the `--rebuild` option');
      rebuild = true;
    }
  }

  color.log('Starting package checks...');
  const checkStart = Date.now();
  const jsTsFiles = files.filter((f) => !f.startsWith('dist/') && (f.endsWith('.js') || f.endsWith('.ts')));
  const packageFilePath = path.join(curDir, 'package.json');
  const json = JSON.parse(fs.readFileSync(packageFilePath, {encoding: 'utf-8'}));

  if (isWebpack) {
    const webpackConfigPath = path.join(curDir, 'webpack.config.js');
    const content = fs.readFileSync(webpackConfigPath, {encoding: 'utf-8'});
  }

  const funcFiles = jsTsFiles.filter((f) => packageFiles.includes(f));
  color.log(`Checks finished in ${Date.now() - checkStart} ms`);
  const reg = new RegExp(/\${(\w*)}/g);
  const errs: string[] = [];

  files.forEach((file: any) => {
    const fullPath = file;
    const relativePath = path.relative(curDir, fullPath);
    const canonicalRelativePath = relativePath.replace(/\\/g, '/');
    if (canonicalRelativePath.includes('/.'))
      return;
    if (canonicalRelativePath.startsWith('.'))
      return;
    if (relativePath.startsWith('node_modules'))
      return;
    if (relativePath.startsWith('dist') && rebuild)
      return;
    if (relativePath.startsWith('src') && !rebuild && isWebpack) {
      if (!relativePath.startsWith('src/package') && !relativePath.startsWith('src\\package'))
        return;
    }
    if (relativePath.startsWith('upload.keys.json'))
      return;
    if (relativePath === 'zip')
      return;
    const t = fs.statSync(fullPath).mtime.toUTCString();
    localTimestamps[canonicalRelativePath] = t;
    if (debug && timestamps[canonicalRelativePath] === t) {
      color.log(`Skipping ${canonicalRelativePath}`);
      return;
    }
    if (canonicalRelativePath.startsWith('connections/')) {
      let f = fs.readFileSync(fullPath, {encoding: 'utf-8'});
      const matches = [...f.matchAll(reg)];
      for (const m of matches) {
        const envVar = process.env[m[1]];
        if (!envVar) {
          errs.push(`${canonicalRelativePath}: cannot find environment variable "${m[1]}"`);
          continue;
        }
        f = f.replace(m[0], envVar);
      }
      zip.append(f, {name: relativePath});
      color.log(`Adding ${canonicalRelativePath}...`);
      return;
    }
    zip.append(fs.createReadStream(fullPath), {name: relativePath});
    color.log(`Adding ${canonicalRelativePath}...`);
  });
  if (errs.length) {
    errs.forEach((e) => color.error(e));
    return 1;
  }

  zip.append(JSON.stringify(localTimestamps), {name: 'timestamps.json'});

  // Upload
  url += `?debug=${debug.toString()}&rebuild=${rebuild.toString()}&dropDb=${(dropDb ?? false).toString()}`;
  if (suffix)
    url += `&suffix=${suffix.toString()}`;
  await zip.finalize();
  const zipBuffer = Buffer.concat(chunks);

  try {
    const body = await fetch(url, {
      method: 'POST',
      body: zipBuffer,
    });
    const log = JSON.parse(await body.text());

    fs.unlinkSync('zip');
    if (log != undefined) {
      if (log['#type'] === 'ApiError') {
        color.error(log['message']);
        console.error(log['innerMessage']);
        console.log(log);
        return 1;
      } else {
        if (color.isVerbose())
          console.log(log);
        color.success(`✓ Published to ${hostAlias || new URL(host).hostname}`);
      }
    }
  } catch (error) {
    if (utils.isConnectivityError(error))
      color.error(`Server is possibly offline: ${url}`);
    if (color.isVerbose())
      console.error(error);
    return 1;
  }
  return 0;
}

export async function publish(args: PublishArgs) {
  color.setVerbose(args.verbose || args.v || false);
  if (!args['skip-check'] && !args['all']&& !args['refresh'])
    check({_: ['check'], verbose: args.verbose, v: args.v});
  config = yaml.load(fs.readFileSync(confPath, {encoding: 'utf-8'})) as utils.Config;
  if (args.refresh) {
    if (path.basename(curDir) !== 'packages')
      curDir = path.dirname(curDir);
    let host = config.default;
    if (args['_'].length === 2)
      host = args['_'][1];
    utils.setHost(host, config);

    const baseUrl = config['servers'][host]['url'];
    let url: string = process.env.HOST ?? '';
    const cfg = getDevKey(url);
    url = cfg.url;

    const key = cfg.key;
    const token = await getToken(url, key);

    url = `${baseUrl}/packages/published/current`;

    let packagesToLoad = ['all'];
    if (args.refresh) {
      packagesToLoad = (await (await fetch(url, {
        method: 'GET',
        headers: {
          'Authorization': token,
        },
      })).json()).map((item: any) => item.name);
    }
    color.log('Loading packages:');
    await loadPackages(curDir, packagesToLoad.join(' '), host, false, false, args.link, args.release);
  } else {
    if (args.link) {
      color.log('Linking');

      await utils.runScript(`npm install`, curDir);
      await utils.runScript(`grok link`, curDir);
      await utils.runScript(`npm run build`, curDir);
    }
    await publishPackage(args);
  }
}

async function publishPackage(args: PublishArgs) {
  const nArgs = args['_'].length;

  if (!args.link) {
    if (args.build || args.rebuild) {
      color.log('Building');
      await utils.runScript('npm install', curDir, false);
      await utils.runScript('npm run build', curDir, false);
    }
  }

  if (args.debug && args.release) {
    color.error('Incompatible options: --debug and --release');
    return false;
  }

  color.log('Publishing...');
  // Create `config.yaml` if it doesn't exist yet
  if (!fs.existsSync(grokDir)) fs.mkdirSync(grokDir);
  if (!fs.existsSync(confPath)) fs.writeFileSync(confPath, yaml.dump(confTemplate));

  let host = config.default;
  const urls = utils.mapURL(config);
  if (nArgs === 2) host = args['_'][1];
  let key = '';
  let url = '';

  // The host can be passed either as a URL or an alias
  try {
    url = new URL(host).href;
    if (url.endsWith('/')) url = url.slice(0, -1);
    if (url in urls) key = config['servers'][urls[url]]['key'];
  } catch (error) {
    if (!(host in config.servers)) return color.error(`Unknown server alias. Please add it to ${confPath}`);
    url = config['servers'][host]['url'];
    key = config['servers'][host]['key'];
  }

  // Update the developer key
  if (args.key) key = args.key;
  if (key === '') return color.warn('Please provide the key with `--key` option or add it by running `grok config`');

  // Get the package name
  if (!fs.existsSync(packDir)) return color.error('`package.json` doesn\'t exist');
  const _package = JSON.parse(fs.readFileSync(packDir, {encoding: 'utf-8'}));
  const packageName = _package.name;

  // Upload the package
  process.env.NODE_TLS_REJECT_UNAUTHORIZED = '0';
  process.on('beforeExit', async () => {
    let code = 0;
    try {

      exec('git rev-parse HEAD', (error: any, stdout: any, stderr: any) => {
        if (error) {
          console.error(`Error executing command: ${error.message}`);
          return;
        }
        if (stderr) {
          console.error(`Standard Error: ${stderr}`);
          return;
        }
        if (!args.suffix && stdout)
          args.suffix = stdout.toString().substring(0, 8);
      });
      await utils.delay(100);
      code = await processPackage(!args.release, Boolean(args.rebuild), url, key, packageName, args.dropDb ?? false, args.suffix, host);
    } catch (error) {
      console.error(error);
      code = 1;
    }
    color.log(`Exiting with code ${code}`);
    process.exit(code);
  });

  return true;
}

interface PublishArgs {
  _: string[],
  build?: boolean,
  ['skip-check']?: boolean,
  rebuild?: boolean,
  debug?: boolean,
  release?: boolean,
  key?: string,
  suffix?: string,
  all?: boolean,
  refresh?: boolean,
  link?: boolean,
  dropDb?: boolean,
  verbose?: boolean,
  v?: boolean
}
