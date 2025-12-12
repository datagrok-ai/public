import fs from 'fs';
import path from 'path';
import walk from 'ignore-walk';
import * as utils from '../utils/utils';
import * as color from '../utils/color-utils';
import {FuncMetadata, FuncParam, FuncValidator, ValidationResult} from '../utils/interfaces';
import {PackageFile} from '../utils/interfaces';
import * as testUtils from '../utils/test-utils';
import {error} from 'console';
import {FuncRoleDescription, functionRoles} from 'datagrok-api/src/const';


const warns = ['Latest package version', 'Datagrok API version should contain'];
const forbiddenNames = ['function', 'class', 'export'];
const namesInFiles = new Map<string, string[]>();

export function check(args: CheckArgs): boolean {
  const nOptions = Object.keys(args).length - 1;
  const curDir = process.cwd();

  if (args.recursive)
    return runChecksRec(curDir, args.soft ?? false);
  else {
    if (!utils.isPackageDir(curDir)) {
      color.error('File `package.json` not found. Run the command from the package directory');
      return false;
    }
    return runChecks(curDir, args.soft ?? false);
  }
}

function runChecks(packagePath: string, soft: boolean = false): boolean {
  if (packagePath.includes(`${path.sep}node_modules${path.sep}`))
    return true;
  const files = (walk.sync({path: packagePath, ignoreFiles: ['.npmignore', '.gitignore']})).filter((e) => !e.includes('node_modules'));
  const jsTsFiles = files.filter((f) => (f.startsWith('src' + path.sep) || f.startsWith('queries' + path.sep) || f.startsWith('scripts' + path.sep)) && (f.endsWith('.js') || f.endsWith('.ts') || f.endsWith('.sql') || f.endsWith('.py')));
  const packageFiles = ['src/package.ts', 'src/detectors.ts', 'src/package.js', 'src/detectors.js',
    'src/package-test.ts', 'src/package-test.js', 'package.js', 'detectors.js'];
  // const funcFiles = jsTsFiles.filter((f) => packageFiles.includes(f)); 
  const errors: string[] = [];
  const warnings: string[] = [];
  const packageFilePath = path.join(packagePath, 'package.json');
  const json: PackageFile = JSON.parse(fs.readFileSync(packageFilePath, {encoding: 'utf-8'}));

  const webpackConfigPath = path.join(packagePath, 'webpack.config.js');
  const isWebpack = fs.existsSync(webpackConfigPath);
  let isReleaseCandidateVersion: boolean = false;
  let externals: { [key: string]: string } | null = null;

  if (/\d+.\d+.\d+-rc(.[A-Za-z0-9]*.[A-Za-z0-9]*)?/.test(json.version))
    isReleaseCandidateVersion = true;
  if (isWebpack) {
    const content = fs.readFileSync(webpackConfigPath, {encoding: 'utf-8'});
    externals = extractExternals(content);
    if (externals)
      errors.push(...checkImportStatements(packagePath, jsTsFiles, externals));
  }
  if (!soft)
    errors.push(...checkSourceMap(packagePath));
  errors.push(...checkNpmIgnore(packagePath));
  warnings.push(...checkScriptNames(packagePath));
  const [signatureWarnings, signatureErrors] = checkFuncSignatures(packagePath, jsTsFiles);
  warnings.push(...signatureWarnings);
  errors.push(...signatureErrors);
  errors.push(...checkPackageFile(packagePath, json, {isWebpack, externals, isReleaseCandidateVersion}));
  if (!isReleaseCandidateVersion)
    warnings.push(...checkChangelog(packagePath, json));

  if (warnings.length) {
    console.log(`${path.basename(packagePath)} warnings`);
    warn(warnings);
  }

  if (errors.length) {
    console.log(`Checking package ${path.basename(packagePath)}...`);
    showError(errors);
    if (soft || json.version.startsWith('0') || (errors.every((w) => warns.some((ww) => w.includes(ww)))))
      return true;
    testUtils.exitWithCode(1);
  }
  console.log(`Checking package ${path.basename(packagePath)}...\t\t\t\u2713 OK`);
  return true;
}

function runChecksRec(dir: string, soft: boolean = false): boolean {
  const files = fs.readdirSync(dir);
  for (const file of files) {
    const filepath = path.join(dir, file);
    const stats = fs.statSync(filepath);
    if (stats.isDirectory()) {
      if (utils.isPackageDir(filepath))
        return runChecks(filepath, soft);
      else {
        if (file !== 'node_modules' && !file.startsWith('.'))
          runChecksRec(path.join(dir, file), soft);
      }
    }
  }
  return false;
}

export function extractExternals(config: string): {} | null {
  const externalsRegex = /(?<=externals)\s*:\s*(\{[\S\s]*?\})/;
  const match = config.match(externalsRegex);
  if (match) {
    // Replace single quotes, comments, and a trailing comma to make a string JSON-like
    const externalStr = match[1]
      .replace(/'/g, '"')
      .replace(/\/\/.*(\r\n|\r|\n)/, '')
      .replace(/(?<=[\S\s]),(?=\s*\})/, '');
    try {
      const externals = JSON.parse(externalStr);
      return externals;
    } catch (e) {
      return null;
    }
  }
  return null;
}

export function checkImportStatements(packagePath: string, files: string[], externals: {}): string[] {
  const modules = [];
  for (const key in externals) {
    modules.push(key);
    if (key.includes('/'))
      modules.push(key.split('/', 1)[0]);
  }
  const importRegex = new RegExp(`^(?!\\/{2})\\s*import\\s+.*(${modules.join('|')}).*(?=\\s+?)`, 'g');
  const validImportRegex = new RegExp(`import\\s+.*(${Object.keys(externals).join('|')})['"]{1}`);
  const warnings: string[] = [];

  function validateImport(file: string, s: string): ValidationResult {
    const value = validImportRegex.test(s);
    const message = value ? '' : 'Pay attention to file ' + file + ': import statement `' +
      s + '` differs from the path given in the webpack config as an external module. ' +
      'It can increase the bundle size.';
    return {value, message};
  }

  for (const file of files) {
    const content = fs.readFileSync(path.join(packagePath, file), {encoding: 'utf-8'});
    const matchedImports = content.match(importRegex);
    if (matchedImports) {
      for (const match of matchedImports) {
        const vr = validateImport(file, match);
        if (!vr.value)
          warnings.push(vr.message);
      }
    }
  }

  return warnings;
}

function normalizeType(type: string): string {
  return type.toLowerCase().replace(/_/g, '');
}

function typesMatch(actual: string, expected: string): boolean {
  return normalizeType(actual) === normalizeType(expected);
}

function parseSignature(sig: string): { inputs: FuncParam[]; outputs: FuncParam[]; variadicIndex: number } {
  const match = sig.match(/^[^(]+\(([^)]*)\)\s*:\s*(.+)$/);
  if (!match) throw new Error(`Invalid signature format: ${sig}`);

  const paramsStr = match[1].trim();
  const returnTypeStr = match[2].trim();

  let variadicIndex = -1;

  const parseParam = (p: string, index: number): FuncParam => {
    const trimmed = p.trim();
    const isVariadic = trimmed.startsWith('...');
    if (isVariadic) variadicIndex = index;

    const paramBody = isVariadic ? trimmed.slice(3) : trimmed;
    const [name, type] = paramBody.split(':').map((s) => s.trim());
    return {name, type: type || 'any'};
  };

  const inputs: FuncParam[] = paramsStr ? paramsStr.split(',').map((p, i) => parseParam(p, i)) : [];
  const outputs: FuncParam[] = returnTypeStr ? returnTypeStr.split('|').map((t) => ({name: '', type: t.trim()})) : [];

  return {inputs, outputs, variadicIndex};
}

function validateFunctionSignature(func: FuncMetadata, roleDesc: FuncRoleDescription): ValidationResult {
  let valid = true;
  let message = '';

  const addError = (msg: string) => {
    valid = false;
    message += msg + '\n';
  };

  if (roleDesc.role === 'app' && func.name) {
    const lower = func.name.toLowerCase();
    if (lower.startsWith('app'))
      addError('Prefix "App" is not needed. Consider removing it.');
    if (lower.endsWith('app'))
      addError('Postfix "App" is not needed. Consider removing it.');
  }

  if (roleDesc.role === 'fileExporter' && (!func.description || func.description === '')) 
    addError('File exporters should have a description parameter');
  
  if (roleDesc.role === 'fileViewer') {
    if (!func.tags || func.tags.length !== 1 || func.tags[0] !== 'fileViewer') 
      addError('File viewers must have only one tag: "fileViewer"');
  }

  const parsed = parseSignature(roleDesc.signature!);
  const maxIndex = parsed.variadicIndex >= 0 ? parsed.variadicIndex : parsed.inputs.length;

  for (let i = 0; i < maxIndex; i++) {
    const expected = parsed.inputs[i];
    const actual = func.inputs[i];

    if (!actual) {
      addError(`Input ${i} missing`);
      continue;
    }

    if (expected.type!.toLowerCase() !== 'any' && !typesMatch(actual.type!, expected.type!)) 
      addError(`Input ${i} type mismatch: expected ${expected.type}, got ${actual.type}`);
  }

  parsed.outputs.forEach((expected, i) => {
    const actual = func.outputs[i];
    if (!actual) 
      addError(`Output ${i} missing`);
    else if (expected.type!.toLowerCase() !== 'any' && !typesMatch(actual.type!, expected.type!)) 
      addError(`Output ${i} type mismatch: expected ${expected.type}, got ${actual.type}`);
    
  });

  [...func.inputs, ...func.outputs].forEach((p) => {
    if (!p.name || !p.type) 
      addError(`Parameter missing name or type: ${p.name ?? '<unnamed>'} (${p.type ?? '<undefined>'})`);
    
  });

  return {value: valid, message};
}

export function checkFuncSignatures(packagePath: string, files: string[]): [string[], string[]] {
  const warnings: string[] = [];
  const errors: string[] = [];
  const namesInFiles = new Map<string, string[]>();
  const roleMap = new Map(functionRoles.map((r) => [r.role, r]));

  for (const file of files) {
    if (file.includes('.min.')) continue;

    const content = fs.readFileSync(path.join(packagePath, file), 'utf-8');
    const functions = getFuncMetadata(content, file.split('.').pop() ?? 'ts');

    for (const f of functions.meta) {
      const allRoles = [
        ...(f.tags || []),
        ...(f.meta?.roles || []),
      ];

      const roles = Array.from(roleMap.keys()).filter((role) => allRoles.includes(role));

      if (roles.length > 1) warnings.push(`File ${file}, function ${f.name}: several function roles are used (${roles.join(', ')})`);
      else if (roles.length === 1) {
        const roleDesc = roleMap.get(roles[0])!;
        const vr = validateFunctionSignature(f, roleDesc);
        if (!vr.value) warnings.push(`File ${file}, function ${f.name}:\n${vr.message}`);
      }

      const invalidNames = f.inputs.filter((e) => forbiddenNames.includes(e?.name ?? ''));
      if (invalidNames.length)
        errors.push(`File ${file}, function ${f.name}: Wrong input names: (${invalidNames.map((e) => e.name).join(', ')})`);

      if (f.name && f.name !== 'postprocess') {
        if (!namesInFiles.has(f.name)) namesInFiles.set(f.name, []);
        namesInFiles.get(f.name)!.push(file);
      }

      if (f.isInvalidateOnWithoutCache)
        errors.push(`File ${file}, function ${f.name}: Can't use invalidateOn without cache`);
      if (f.cache && !utils.cacheValues.includes(f.cache))
        errors.push(`File ${file}, function ${f.name}: unsupported cache variable: ${f.cache}`);
      if (f.invalidateOn && !utils.isValidCron(f.invalidateOn))
        errors.push(`File ${file}, function ${f.name}: unsupported invalidateOn: ${f.invalidateOn}`);
    }

    functions.warnings.forEach((w) => warnings.push(`${w} In the file: ${file}.`));
  }

  for (const [name, files] of namesInFiles) {
    if (files.length > 1)
      errors.push(`Duplicate names ('${name}'): \n  ${files.join('\n  ')}`);
  }

  return [warnings, errors];
}

const sharedLibExternals: { [lib: string]: {} } = {
  'common/html2canvas.min.js': {'exceljs': 'ExcelJS'},
  'common/exceljs.min.js': {'html2canvas': 'html2canvas'},
  'common/ngl_viewer/ngl.js': {'NGL': 'NGL'},
  'common/openchemlib-full.js': {'openchemlib/full': 'OCL'},
  'common/codemirror/codemirror.js': {'codemirror': 'CodeMirror'},
  'common/vue.js': {'vue': 'Vue'},
};

export function checkPackageFile(packagePath: string, json: PackageFile, options?: {
  externals?:
  { [key: string]: string } | null, isWebpack?: boolean | null, isReleaseCandidateVersion?: boolean
}): string[] {
  const warnings: string[] = [];
  const isPublicPackage = path.basename(path.dirname(packagePath)) === 'packages' &&
    path.basename(path.dirname(path.dirname(packagePath))) === 'public';

  if (!json.description)
    warnings.push('File "package.json": "description" field is empty. Provide a package description.');

  if (Array.isArray(json.properties) && json.properties.length > 0) {
    for (const propInfo of json.properties) {
      if (typeof propInfo !== 'object')
        warnings.push('File "package.json": Invalid property annotation in the "properties" field.');
      else if (!propInfo.name)
        warnings.push('File "package.json": Add a property name for each property in the "properties" field.');
      else if (!utils.propertyTypes.includes(propInfo.propertyType))
        warnings.push(`File "package.json": Invalid property type for property ${propInfo.name}.`);
    }
  }

  if (json.repository == null && isPublicPackage)
    warnings.push('File "package.json": add the "repository" field.');

  if (json.author == null && isPublicPackage)
    warnings.push('File "package.json": add the "author" field.');

  if (json.version.includes('beta') && isPublicPackage)
    warnings.push('File "package.json": public package cannot be beta version.');

  const api = json.dependencies?.['datagrok-api'];
  if (api) {
    if (api === '../../js-api') { } else if (api === 'latest')
      warnings.push('File "package.json": you should specify Datagrok API version constraint (for example ^1.16.0, >=1.16.0).');
    else if (options?.isReleaseCandidateVersion === false && (!/^(\^|>|<|~).+/.test(api)))
      warnings.push('File "package.json": Datagrok API version should starts with > | >= | ~ | ^ | < | <=');
  }

  // const dt = json.devDependencies?.['datagrok-tools'] ?? json.dependencies?.['datagrok-tools'];
  // if (dt && dt !== 'latest')
  //   warnings.push('File "package.json": "datagrok-tools" dependency must be "latest" version.');

  if (Array.isArray(json.sources) && json.sources.length > 0) {
    for (const source of json.sources) {
      if (typeof source !== 'string')
        warnings.push(`File "package.json": Only file paths and URLs are allowed in sources. Modify the source ${source}`);
      if (utils.absUrlRegex.test(source))
        continue;
      if (source.startsWith('common/')) {
        if (options?.isWebpack && source.endsWith('.js')) {
          if (options?.externals) {
            if (source in sharedLibExternals) {
              const [lib, name] = Object.entries(sharedLibExternals[source])[0];
              if (!(lib in options.externals && options.externals[lib] === name)) {
                warnings.push(`Webpack config parsing: Consider adding source "${source}" to webpack externals:\n` +
                  `'${lib}': '${name}'\n`);
              }
            } else
              warnings.push(`File "package.json": source "${source}" not in the list of shared libraries`);

          } else {
            warnings.push('Webpack config parsing: External modules not found.\n' +
              `Consider adding source "${source}" to webpack externals` + (source in sharedLibExternals ? ':\n' +
                `'${Object.keys(sharedLibExternals[source])[0]}': '${Object.values(sharedLibExternals[source])[0]}'\n` : ''));
          }
        }
        continue;
      }
      if (source.startsWith('src/') && fs.existsSync(path.join(packagePath, 'webpack.config.js'))) {
        warnings.push('File "package.json": Sources cannot include files from the \`src/\` directory. ' +
          `Move file ${source} to another folder.`);
      }
      if (!(fs.existsSync(path.join(packagePath, source))))
        warnings.push(`Source ${source} not found in the package.`);
    }
  }

  if (options?.isReleaseCandidateVersion === true) {
    let hasRCDependency = false;

    for (const dependency of Object.keys(json.dependencies ?? {})) {
      if (/\d+.\d+.\d+-rc(.[A-Za-z0-9]*.[A-Za-z0-9]*)?/.test((json.dependencies ?? {})[dependency])) {
        hasRCDependency = true;
        break;
      }
    }

    if (!hasRCDependency) {
      for (const dependency of Object.keys(json.dependencies ?? {})) {
        if (/\d+.\d+.\d+-rc(.[A-Za-z0-9]*.[A-Za-z0-9]*)?/.test((json.devDependencies ?? {})[dependency])) {
          hasRCDependency = true;
          break;
        }
      }
    }

    if (!hasRCDependency)
      warnings.push('Release candidate error: Current package doesnt have any dependencies from any release candidate package ');
  }

  return warnings;
}

export function checkChangelog(packagePath: string, json: PackageFile): string[] {
  if (json.servicePackage) return [];
  const warnings: string[] = [];
  let clf: string;
  try {
    clf = fs.readFileSync(path.join(packagePath, 'CHANGELOG.md'), {encoding: 'utf-8'});
  } catch (e) {
    return ['CHANGELOG.md file does not exist\n'];
  }
  let regex = /^##[^#].*$/gm;
  const h2 = clf.match(regex);
  if (!h2) return ['No versions found in CHANGELOG.md\n'];
  regex = /^## \d+\.\d+\.\d+ \((\d{4}-\d{2}-\d{2}|WIP)\)$/;
  for (const h of h2) {
    if (!regex.test(h))
      warnings.push(`CHANGELOG: '${h}' does not match the h2 format, expected: ## <version> (<yyyy-mm-dd> | WIP)\n`);
  }
  regex = /^## (\d+\.\d+\.\d+)/;
  const v1 = h2[0].match(regex)?.[1];
  const v2 = h2[1]?.match(regex)?.[1];
  if (v1 !== json.version && v2 !== json.version)
    warnings.push(`Latest package version (${json.version}) is not in CHANGELOG\n`);

  if (warnings.length)
    warnings.push('Changelog guideline: https://datagrok.ai/help/develop/dev-process/changelog-policy#changelog-guideline');

  return warnings;
}

export function checkSourceMap(packagePath: string): string[] {
  const warnings: string[] = [];
  const tsconfigFilePath = path.join(packagePath, 'tsconfig.json');
  const webpackConfigFilePath = path.join(packagePath, 'webpack.config.js');

  if (fs.existsSync(tsconfigFilePath) && fs.existsSync(webpackConfigFilePath)) {
    const configJson: string = fs.readFileSync(tsconfigFilePath, {encoding: 'utf-8'}); // cant convert to json because file contains comments

    if (!(new RegExp('"sourceMap"\\s*:\\s*true')).test(configJson))
      warnings.push('ts config doesnt contain source map');

    const webpackConfigJson: string = fs.readFileSync(webpackConfigFilePath, {encoding: 'utf-8'}); // cant convert to json because file contains comments

    if (!(new RegExp(`devtool\\s*:\\s*(([^\\n]*?[^\\n]*source-map[^\\n]*:[^\\n]*source-map[^\\n]*)|('(inline-)?source-map'))\\s*`)).test(webpackConfigJson))
      warnings.push('webpack config doesnt contain source map');

    if (!fs.existsSync(packagePath + '/dist/package.js'))
      warnings.push('dist\\package.js file doesnt exists');

    if (!fs.existsSync(packagePath + '/dist/package-test.js'))
      warnings.push('dist\\package-test.js file doesnt exists');

  }
  return warnings;
}


export function checkNpmIgnore(packagePath: string): string[] {
  const warnings: string[] = [];
  if (fs.existsSync(path.join(...[packagePath, '.npmignore']))) {
    const npmIgnoreContent: string = fs.readFileSync(path.join(...[packagePath, '.npmignore']), {encoding: 'utf-8'});
    for (const row of npmIgnoreContent.split('\n')) {
      if ((row.match(new RegExp('\\s*dist\\/?\\s*$'))?.length || -1) > 0) {
        warnings.push('there is dist directory in .npmignore');
        break;
      }
    }
  } else
    warnings.push('.npmignore doesnt exists');

  return warnings;
}

function checkScriptNames(packagePath: string): string[] {
  const warnings: string[] = [];

  try {
    if (fs.existsSync(packagePath)) {
      const filesInDirectory = getAllFilesInDirectory(packagePath);
      for (const fileName of filesInDirectory) {
        if (fileName.match(new RegExp('^[A-Za-z0-9._-]*$'))?.length !== 1)
          warnings.push(`${fileName}: file name contains inappropriate symbols`);
      }
    }
  } catch (err) {
  }
  return warnings;
}

function getAllFilesInDirectory(directoryPath: string): string[] {
  const excludedFilesToCheck: string[] = ['node_modules', 'dist'];

  let fileNames: string[] = [];
  const entries = fs.readdirSync(directoryPath);
  entries.forEach((entry) => {
    const entryPath = path.join(directoryPath, entry);
    const stat = fs.statSync(entryPath);

    if (stat.isFile()) 
      fileNames.push(entry);
    else if (stat.isDirectory() && !excludedFilesToCheck.includes(entry)) {
      const subDirectoryFiles = getAllFilesInDirectory(entryPath);
      fileNames = fileNames.concat(subDirectoryFiles);
    }
  });
  return fileNames;
}

function showError(errors: string[]): void {
  errors.forEach((w) => color.error(w));
}

function warn(warnings: string[]): void {
  warnings.forEach((w) => color.warn(w));
}

function getFuncMetadata(script: string, fileExtention: string): { meta: FuncMetadata[], warnings: string[] } {
  const funcData: FuncMetadata[] = [];
  const warnings: string[] = [];
  let isHeader = false;
  let data: FuncMetadata = {name: '', inputs: [], outputs: []};

  for (const line of script.split('\n')) {
    if (!line)
      continue;
    //@ts-ignore
    const match = line.match(utils.fileParamRegex[fileExtention]);
    if (match) {
      if (!isHeader)
        isHeader = true;
      const param = match[1];

      if (!utils.headerTags.includes(param) && !param.includes('meta.')) {
        warnings.push(`Unknown header tag: ${param},`);
        continue;
      }
      if (param === 'name')
        data.name = line.match(utils.nameAnnRegex)?.[2]?.toLocaleLowerCase();
      else if (param === 'description')
        data.description = match[2];
      else if (param === 'input') 
        data.inputs.push({type: match[2], name: match[3]});
      
      else if (param === 'output')
        data.outputs.push({type: match[2], name: match[3]});
      else if (param === 'tags') 
        data.tags = match.input && match[3] ? match.input.split(':')[1].split(',').map((t) => t.trim()) : [match[2]];

      else if (param === 'meta.roles') {
        data.meta = data.meta || {};
        data.meta.roles = match[2].split(',').map((role) => role.trim());
      } else if (param === 'meta.cache') 
        data.cache = line.split(':').pop()?.trim();
      
      else if (param === 'meta.cache.invalidateOn') 
        data.invalidateOn = line.split(':').pop()?.trim();
      
      else if (param === 'meta.invalidateOn') 
        data.isInvalidateOnWithoutCache = true;
      
    }
    if (isHeader) {
      const nm = line.match(utils.nameRegex);
      if (data.name === '') {
        data = {name: '', inputs: [], outputs: []};
        isHeader = false;
        continue;
      }
      if (nm)
        data.name = nm[1]?.toLocaleLowerCase();
      if (data.name && !match) {
        funcData.push(data);
        data = {name: '', inputs: [], outputs: []};
        isHeader = false;
      }
    }
  }

  return {meta: funcData, warnings};
}

interface CheckArgs {
  _: string[],
  r?: boolean,
  recursive?: boolean,
  soft?: boolean,
}
