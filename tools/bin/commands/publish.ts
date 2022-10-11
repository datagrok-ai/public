// @ts-ignore
import archiver from 'archiver-promise';
import fs from 'fs';
// @ts-ignore
import fetch from 'node-fetch';
import os from 'os';
import path from 'path';
import walk from 'ignore-walk';
import yaml from 'js-yaml';
import { checkImportStatements, checkFuncSignatures, extractExternals } from './check';
import * as utils from '../utils/utils';
import { Indexable } from '../utils/utils';
import * as color from '../utils/color-utils';


const grokDir = path.join(os.homedir(), '.grok');
const confPath = path.join(grokDir, 'config.yaml');
const confTemplateDir = path.join(path.dirname(path.dirname(__dirname)), 'config-template.yaml');
const confTemplate = yaml.load(fs.readFileSync(confTemplateDir, { encoding: 'utf-8' }));
const curDir = process.cwd();
const packDir = path.join(curDir, 'package.json');
const packageFiles = [
  'src/package.ts', 'src/detectors.ts', 'src/package.js', 'src/detectors.js',
  'src/package-test.ts', 'src/package-test.js', 'package.js', 'detectors.js',
];

export async function processPackage(debug: boolean, rebuild: boolean, host: string, devKey: string, packageName: any, suffix?: string) {
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
      console.error(error);
      return 1;
    }
  }

  let zip = archiver('zip', {store: false});

  // Gather the files
  let localTimestamps: Indexable = {};
  let files = await walk({
    path: '.',
    ignoreFiles: ['.npmignore', '.gitignore', '.grokignore'],
    includeEmpty: false,
    follow: true
  });

  let isWebpack = fs.existsSync('webpack.config.js');
  if (!rebuild && isWebpack) {
    if (fs.existsSync('dist/package.js')) {
      const distFiles = await walk({
        path: './dist',
        ignoreFiles: [],
        includeEmpty: false,
        follow: true
      });
      distFiles.forEach((df: any) => {
        files.push(`dist/${df}`);
      })
    } else {
      color.warn('File `dist/package.js` not found. Building the package on the server side...\n' +
        'Next time, please build your package locally with Webpack beforehand\n' +
        'or run `grok publish` with the `--rebuild` option');
      rebuild = true;
    }
  }

  let contentValidationLog = '';
  console.log('Starting package checks...');
  const checkStart = Date.now();
  const jsTsFiles = files.filter((f) => !f.startsWith('dist/') && (f.endsWith('.js') || f.endsWith('.ts')));

  if (isWebpack) {
    const webpackConfigPath = path.join(curDir, 'webpack.config.js');
    const content = fs.readFileSync(webpackConfigPath, { encoding: 'utf-8' });
    const externals = extractExternals(content);
    if (externals) {
      const importWarnings = checkImportStatements(curDir, jsTsFiles, externals);
      contentValidationLog += importWarnings.join('\n') + (importWarnings.length ? '\n' : '');
    }
  }

  const funcFiles = jsTsFiles.filter((f) => packageFiles.includes(f));
  const funcWarnings = checkFuncSignatures(curDir, funcFiles);
  contentValidationLog += funcWarnings.join('\n') + (funcWarnings.length ? '\n' : '');
  console.log(`Checks finished in ${Date.now() - checkStart} ms`);

  files.forEach((file: any) => {
    let fullPath = file;
    let relativePath = path.relative(curDir, fullPath);
    let canonicalRelativePath = relativePath.replace(/\\/g, '/');
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
    if (!utils.checkScriptLocation(canonicalRelativePath)) {
      contentValidationLog += `Warning: file \`${canonicalRelativePath}\`` +
        ` should be in directory \`${path.basename(curDir)}/scripts/\`\n`;
    }
    let t = fs.statSync(fullPath).mtime.toUTCString();
    localTimestamps[canonicalRelativePath] = t;
    if (debug && timestamps[canonicalRelativePath] === t) {
      console.log(`Skipping ${canonicalRelativePath}`);
      return;
    }
    zip.append(fs.createReadStream(fullPath), {name: relativePath});
    console.log(`Adding ${canonicalRelativePath}...`);
  });
  zip.append(JSON.stringify(localTimestamps), {name: 'timestamps.json'});

  // Upload
  url += `?debug=${debug.toString()}&rebuild=${rebuild.toString()}`;
  if (suffix) url += `&suffix=${suffix.toString()}`;
  let uploadPromise = new Promise((resolve, reject) => {
    fetch(url, {
      method: 'POST',
      body: zip
    }).then(async (body: any) => {
      let response;
      try {
        response = await body.text();
        return JSON.parse(response);
      } catch (error) {
        console.error(response);
      }
    }).then((j: {}) => resolve(j)).catch((err: Error) => {
      reject(err);
    });
  }).catch(error => {
    console.error(error);
  });
  await zip.finalize();

  try {
    let log = await uploadPromise as Indexable;

    fs.unlinkSync('zip');
    if (log['#type'] === 'ApiError') {
      color.error(log['message']);
      console.error(log['innerMessage']);
      return 1;
    } else {
      console.log(log);
      color.warn(contentValidationLog);
    }
  } catch (error) {
    console.error(error);
    return 1;
  }
  return 0;
}

export function publish(args: PublishArgs) {
  const nOptions = Object.keys(args).length - 1;
  const nArgs = args['_'].length;

  if (nArgs > 2 || nOptions > 5) return false;
  if (!Object.keys(args).slice(1).every(option => ['build', 'rebuild',
  'debug', 'release', 'k', 'key', 'suffix'].includes(option))) return false;

  if (args.build && args.rebuild) {
    color.error('Incompatible options: --build and --rebuild');
    return false;
  }

  if (args.debug && args.release) {
    color.error('Incompatible options: --debug and --release');
    return false;
  }

  // Create `config.yaml` if it doesn't exist yet
  if (!fs.existsSync(grokDir)) fs.mkdirSync(grokDir);
  if (!fs.existsSync(confPath)) fs.writeFileSync(confPath, yaml.dump(confTemplate));

  let config = yaml.load(fs.readFileSync(confPath, { encoding: 'utf-8' })) as utils.Config;
  let host = config.default;
  let urls = utils.mapURL(config);
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
  let _package = JSON.parse(fs.readFileSync(packDir, { encoding: 'utf-8' }));
  let packageName = _package.name;

  // Upload the package
  process.env.NODE_TLS_REJECT_UNAUTHORIZED = "0";
  process.on('beforeExit', async () => {
    let code = 0;
    try {
      code = await processPackage(!args.release, Boolean(args.rebuild), url, key, packageName, args.suffix)

    } catch (error) {
      console.error(error);
      code = 1;
    }
    console.log(`Exiting with code ${code}`);
    process.exit(code);
  });

  return true;
}

interface PublishArgs {
  _: string[],
  build?: boolean,
  rebuild?: boolean,
  debug?: boolean,
  release?: boolean,
  key?: string,
  suffix?: string,
}
