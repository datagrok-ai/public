// @ts-ignore
import archiver from 'archiver-promise';
import crypto from 'crypto';
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
import {generateCeleryArtifacts} from '../utils/python-celery-gen';

const {exec, execSync} = require('child_process');

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

const BLEEDING_EDGE_TAG = 'bleeding-edge';

interface DockerImageInfo {
  imageName: string;
  imageTag: string;
  dirName: string;
  fullLocalName: string;
}

interface DockerImageResult {
  image: string | null;
  fallback: boolean;
  requestedVersion?: string;
}

function discoverDockerfiles(packageName: string, version: string, debug: boolean): DockerImageInfo[] {
  const dockerfilesDir = path.join(curDir, 'dockerfiles');
  if (!fs.existsSync(dockerfilesDir))
    return [];

  const results: DockerImageInfo[] = [];
  const entries = fs.readdirSync(dockerfilesDir, {withFileTypes: true});
  for (const entry of entries) {
    if (!entry.isDirectory())
      continue;
    const dockerfilePath = path.join(dockerfilesDir, entry.name, 'Dockerfile');
    if (!fs.existsSync(dockerfilePath))
      continue;
    const cleanName = utils.removeScope(packageName).toLowerCase();
    const imageName = `${cleanName}-${entry.name.toLowerCase()}`;
    const imageTag = version;
    results.push({
      imageName,
      imageTag,
      dirName: entry.name,
      fullLocalName: `${imageName}:${imageTag}`,
    });
  }

  // Handle Dockerfile directly in dockerfiles/ (single-container layout)
  const rootDockerfile = path.join(dockerfilesDir, 'Dockerfile');
  if (fs.existsSync(rootDockerfile)) {
    const cleanName = utils.removeScope(packageName).toLowerCase();
    results.push({
      imageName: cleanName,
      imageTag: version,
      dirName: '.',
      fullLocalName: `${cleanName}:${version}`,
    });
  }

  return results;
}

function dockerCommand(args: string): string {
  return execSync(`docker ${args}`, {encoding: 'utf-8', stdio: ['pipe', 'pipe', 'pipe']}).trim();
}

function localImageExists(fullName: string, checkPlatform: boolean = true): boolean {
  try {
    const output = dockerCommand(`images --format "{{.Repository}}:{{.Tag}}"`);
    if (!output.split('\n').some((line: string) => line.trim() === fullName))
      return false;
    if (checkPlatform) {
      // Verify the image is linux/amd64 — ARM images won't run on Datagrok servers
      const inspect = dockerCommand(`inspect --format "{{.Os}}/{{.Architecture}}" ${fullName}`);
      if (inspect.trim() !== 'linux/amd64') {
        color.warn(`Local image ${fullName} is ${inspect.trim()}, not linux/amd64 — treating as not found`);
        return false;
      }
    }
    return true;
  } catch {
    return false;
  }
}

function dockerLogin(registry: string, devKey: string): boolean {
  try {
    dockerCommand(`login ${registry} -u any -p ${devKey}`);
    return true;
  } catch (e: any) {
    color.warn(`Docker login to ${registry} failed: ${e.message || e}`);
    return false;
  }
}

function dockerTag(source: string, target: string): boolean {
  try {
    dockerCommand(`tag ${source} ${target}`);
    return true;
  } catch (e: any) {
    color.error(`Failed to tag ${source} as ${target}: ${e.message || e}`);
    return false;
  }
}

function dockerPush(image: string): boolean {
  try {
    execSync(`docker push ${image}`, {encoding: 'utf-8', stdio: 'inherit'});
    return true;
  } catch (e: any) {
    color.error(`Failed to push ${image}: ${e.message || e}`);
    return false;
  }
}

function isLocalhostRegistry(registry: string): boolean {
  return registry.startsWith('localhost') || registry.startsWith('127.0.0.1');
}

function dockerRemove(imageName: string): boolean {
  try {
    dockerCommand(`rmi ${imageName}`);
    return true;
  }
  catch {
    return false;
  }
}

function dockerBuild(imageName: string, dockerfilePath: string, context: string, usePlatform: boolean = true): boolean {
  try {
    const platformFlag = usePlatform ? '--platform linux/amd64 ' : '';
    color.log(`Building Docker image ${imageName}${usePlatform ? ' (platform: linux/amd64)' : ''}...`);
    execSync(`docker build ${platformFlag}-t ${imageName} -f ${dockerfilePath} ${context}`, {
      encoding: 'utf-8',
      stdio: 'inherit',
    });
    return true;
  } catch (e: any) {
    color.error(`Failed to build ${imageName}: ${e.message || e}`);
    return false;
  }
}

function calculateFolderHash(dirPath: string): string {
  const hash = crypto.createHash('sha256');
  const entries = listRecursive(dirPath, '');
  entries.sort((a, b) => a.relPath < b.relPath ? -1 : a.relPath > b.relPath ? 1 : 0);
  for (const entry of entries) {
    if (entry.isDir) {
      color.log(`  Hash entry: dir:${entry.relPath}:`);
      hash.update(`dir:${entry.relPath}:`);
    }
    else {
      // Normalize CRLF to LF to match server-side storage
      const raw = fs.readFileSync(entry.fullPath);
      const content = raw.includes(0x0d) && !raw.includes(0x00)
        ? Buffer.from(raw.toString('utf8').replace(/\r\n/g, '\n'), 'utf8')
        : raw;
      color.log(`  Hash entry: file:${entry.relPath}: (${content.length} bytes)`);
      hash.update(`file:${entry.relPath}:`);
      hash.update(content);
    }
  }
  const result = hash.digest('hex');
  color.log(`  Folder hash: ${result}`);
  return result;
}

function listRecursive(basePath: string, rel: string): {relPath: string, fullPath: string, isDir: boolean}[] {
  const results: {relPath: string, fullPath: string, isDir: boolean}[] = [];
  const fullDir = rel ? path.join(basePath, rel) : basePath;
  for (const entry of fs.readdirSync(fullDir, {withFileTypes: true})) {
    const relPath = rel ? `${rel}/${entry.name}` : entry.name;
    const fullPath = path.join(fullDir, entry.name);
    if (entry.isDirectory()) {
      results.push({relPath, fullPath, isDir: true});
      results.push(...listRecursive(basePath, relPath));
    }
    else
      results.push({relPath, fullPath, isDir: false});
  }
  return results;
}

async function getUserLogin(host: string, devKey: string): Promise<{login: string, token: string} | null> {
  let loginResp;
  try {
    loginResp = await fetch(`${host}/users/login/dev/${devKey}`, {method: 'POST'});
  } catch (e: any) {
    color.warn(`Cannot reach server ${host}: ${e.message || e}`);
    return null;
  }
  if (loginResp.status !== 200)
    return null;
  const loginData = await loginResp.json();
  const token = loginData.token;
  try {
    const userResp = await fetch(`${host}/users/current`, {headers: {'Authorization': token}});
    if (userResp.status !== 200)
      return null;
    const userData = await userResp.json();
    return {login: userData.login, token};
  } catch (e: any) {
    return null;
  }
}

interface FallbackResult {
  found: {version: string, image: string, hashMatch?: boolean} | null;
  serverError: string | null;
}

async function resolveLatestCompatible(host: string, devKey: string, dockerName: string, version?: string, contentHash?: string): Promise<FallbackResult> {
  const userInfo = await getUserLogin(host, devKey);
  if (!userInfo)
    return {found: null, serverError: `Authentication failed. Check your developer key.`};

  try {
    let url = `${host}/docker/images/${encodeURIComponent(dockerName)}/latest-compatible`;
    const params: string[] = [];
    if (version)
      params.push(`version=${encodeURIComponent(version)}`);
    if (contentHash)
      params.push(`contentHash=${encodeURIComponent(contentHash)}`);
    if (params.length > 0)
      url += '?' + params.join('&');
    const resp = await fetch(url, {headers: {'Authorization': userInfo.token}});
    if (resp.status === 200)
      return {found: await resp.json(), serverError: null};
    if (resp.status === 404)
      return {found: null, serverError: null};
    return {found: null, serverError: `Unexpected response (HTTP ${resp.status}) from latest-compatible endpoint`};
  } catch (e: any) {
    return {found: null, serverError: `Failed to query latest-compatible: ${e.message || e}`};
  }
}

async function processDockerImages(
  packageName: string,
  version: string,
  registry: string | undefined,
  devKey: string,
  host: string,
  rebuildDocker: boolean,
  zip: any,
  localTimestamps: Indexable,
  debug: boolean,
  skipDockerRebuild: boolean = false,
): Promise<void> {
  const dockerImages = discoverDockerfiles(packageName, version, debug);
  if (dockerImages.length === 0)
    return;

  color.log(`Found ${dockerImages.length} Dockerfile(s)`);

  if (registry)
    dockerLogin(registry, devKey);

  const needsPlatform = !!registry && !isLocalhostRegistry(registry);

  for (const img of dockerImages) {
    color.info(`Processing docker image ${img.fullLocalName}...`);
    let result: DockerImageResult;
    const dockerfileDir = path.join('dockerfiles', img.dirName);
    const dockerfilePath = path.join(dockerfileDir, 'Dockerfile');
    const contentHash = calculateFolderHash(path.join(curDir, dockerfileDir));

    // In release mode, append short content hash to tag for cache-busting
    const shortHash = contentHash.substring(0, 8);
    const registryTag = debug ? img.imageTag : `${img.imageTag}.${shortHash}`;
    const remoteFullName = `${img.imageName}:${registryTag}`;

    const buildAndPush = (): DockerImageResult | null => {
      if (dockerBuild(img.fullLocalName, dockerfilePath, dockerfileDir, needsPlatform)) {
        if (registry) {
          const remoteTag = `${registry}/datagrok/${remoteFullName}`;
          dockerTag(img.fullLocalName, remoteTag);
        }
        else
          dockerTag(img.fullLocalName, `datagrok/${remoteFullName}`);
        color.success(`Built and tagged ${img.fullLocalName}`);
        return pushImage(img.imageName, registryTag, registry);
      }
      return null;
    };

    if (rebuildDocker) {
      // Delete old images before rebuilding
      dockerRemove(img.fullLocalName);
      if (registry)
        dockerRemove(`${registry}/datagrok/${remoteFullName}`);
      result = buildAndPush() ?? await fallbackImage(img, host, devKey, registry, version, contentHash);
    }
    else {
      // Look for registry-qualified image first, then unqualified
      let foundLocalName: string | null = null;
      if (registry) {
        const registryQualified = `${registry}/datagrok/${remoteFullName}`;
        if (localImageExists(registryQualified, needsPlatform))
          foundLocalName = registryQualified;
      }
      if (!foundLocalName && localImageExists(img.fullLocalName, needsPlatform))
        foundLocalName = img.fullLocalName;

      if (foundLocalName) {
        color.success(`Found local image ${foundLocalName}`);
        if (registry) {
          const remoteTag = `${registry}/datagrok/${remoteFullName}`;
          if (foundLocalName !== remoteTag)
            dockerTag(foundLocalName, remoteTag);
        }
        else {
          const canonicalTag = `datagrok/${remoteFullName}`;
          dockerTag(foundLocalName, canonicalTag);
          color.log(`  Tagged as ${canonicalTag}`);
        }
        result = pushImage(img.imageName, registryTag, registry);
      }
      else {
        color.warn(`Local image not found. Expected: ${img.fullLocalName}`);
        color.log(`  Build it with: docker build -t ${img.fullLocalName} -f ${dockerfilePath} ${dockerfileDir}`);
        const fallback = await fallbackImage(img, host, devKey, registry, version, contentHash);
        if (fallback.serverError) {
          color.error(`Cannot resolve fallback: ${fallback.serverError}`);
          result = {image: null, fallback: true, requestedVersion: registryTag};
        }
        else if (fallback.image && fallback.hashMatch === true) {
          result = fallback;
          color.success(`Falling back to ${fallback.image} (dockerfile unchanged)`);
        }
        else if (fallback.image && fallback.hashMatch === false && !skipDockerRebuild) {
          color.warn(`Dockerfile folder has changed. Rebuilding image...`);
          result = buildAndPush() ?? {image: fallback.image, fallback: true, requestedVersion: registryTag};
          if (!result || result.fallback)
            color.warn(`Build failed. Falling back to ${fallback.image} (hash mismatch)`);
        }
        else if (skipDockerRebuild) {
          color.warn(`No fallback available. Skipping docker build (--skip-docker-rebuild).`);
          result = {image: null, fallback: true, requestedVersion: registryTag};
        }
        else {
          // No fallback and no local image — must build
          color.warn(`No fallback available. Building ${img.fullLocalName}...`);
          const built = buildAndPush();
          if (built)
            result = built;
          else {
            result = {image: null, fallback: true, requestedVersion: registryTag};
            color.error(`Failed to build ${img.fullLocalName}. No container will be available.`);
          }
        }
      }
    }

    const imageJsonPath = path.join('dockerfiles', img.dirName, 'image.json');
    zip.append(JSON.stringify(result, null, 2), {name: imageJsonPath});
    localTimestamps[imageJsonPath.replace(/\\/g, '/')] = new Date().toUTCString();
    color.log(`Added ${imageJsonPath}`);
  }
}

function pushImage(imageName: string, tag: string, registry: string | undefined): DockerImageResult {
  const canonicalImage = `datagrok/${imageName}:${tag}`;
  if (registry) {
    const remoteTag = `${registry}/${canonicalImage}`;
    // Image should already be tagged from build or retag step
    if (dockerPush(remoteTag))
      return {image: canonicalImage, fallback: false};
    color.warn(`Push failed, image tagged locally only: ${remoteTag}`);
    return {image: canonicalImage, fallback: false};
  }
  color.warn(`No registry configured. Image tagged locally only. Run \`grok config --registry\` to configure.`);
  return {image: canonicalImage, fallback: false};
}

interface FallbackImageResult extends DockerImageResult {
  serverError?: string;
  hashMatch?: boolean;
}

async function fallbackImage(
  img: DockerImageInfo,
  host: string,
  devKey: string,
  registry: string | undefined,
  version?: string,
  contentHash?: string,
): Promise<FallbackImageResult> {
  const {found, serverError} = await resolveLatestCompatible(host, devKey, img.imageName, version, contentHash);
  if (serverError)
    return {image: null, fallback: true, requestedVersion: img.imageTag, serverError};
  if (found)
    return {image: found.image, fallback: true, requestedVersion: img.imageTag, hashMatch: found.hashMatch};
  return {image: null, fallback: true, requestedVersion: img.imageTag};
}

export async function processPackage(debug: boolean, rebuild: boolean, host: string, devKey: string, packageName: any, dropDb: boolean, suffix?: string, hostAlias?: string, registry?: string, rebuildDocker?: boolean, skipDockerRebuild?: boolean) {
  // Validate server connectivity and dev key
  let timestamps: Indexable = {};
  let url = `${host}/packages/dev/${devKey}/${packageName}`;
  try {
    const checkResp = await fetch(url + '/timestamps');
    const checkData = await checkResp.json();
    if (checkData['#type'] === 'ApiError') {
      color.error(checkData.message);
      return 1;
    }
    if (debug)
      timestamps = checkData;
  } catch (error) {
    if (utils.isConnectivityError(error))
      color.error(`Server is possibly offline: ${host}`);
    if (color.isVerbose())
      console.error(error);
    return 1;
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

  // Generate Celery Docker artifacts from python/ if present
  generateCeleryArtifacts(curDir);

  // Process Docker images and inject image.json into the ZIP
  let dockerVersion = json.version;
  if (debug) {
    const userInfo = await getUserLogin(host, devKey);
    if (userInfo)
      dockerVersion = userInfo.login;
  }
  await processDockerImages(packageName, dockerVersion, registry, devKey, host, rebuildDocker ?? false, zip, localTimestamps, debug, skipDockerRebuild ?? false);

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
    check({_: ['check']});
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
    return await publishPackage(args);
  }
  return true;
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
  let registry: string | undefined;

  // The host can be passed either as a URL or an alias
  try {
    url = new URL(host).href;
    if (url.endsWith('/')) url = url.slice(0, -1);
    if (url in urls) {
      const alias = urls[url];
      key = config['servers'][alias]['key'];
      registry = config['servers'][alias]['registry'];
    }
  } catch (error) {
    if (!(host in config.servers)) return color.error(`Unknown server alias. Please add it to ${confPath}`);
    url = config['servers'][host]['url'];
    key = config['servers'][host]['key'];
    registry = config['servers'][host]['registry'];
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
      code = await processPackage(!args.release, Boolean(args.rebuild), url, key, packageName, args.dropDb ?? false, args.suffix, host, registry, args['rebuild-docker'], args['skip-docker-rebuild']);
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
  ['rebuild-docker']?: boolean,
  ['skip-docker-rebuild']?: boolean,
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
