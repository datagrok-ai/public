import fs from 'fs';
import path from 'path';
import walk from 'ignore-walk';
import * as utils from '../utils/utils';
import * as color from '../utils/color-utils';
import { FuncMetadata, FuncParam, FuncValidator, ValidationResult } from '../utils/interfaces';
import { PackageFile } from '../utils/interfaces';
import * as testUtils from '../utils/test-utils';
import { error } from 'console';

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
  const files = walk.sync({ path: packagePath, ignoreFiles: ['.npmignore', '.gitignore'] });
  console.log(files.filter(e=> e.includes('node_modules')))
  const jsTsFiles = files.filter((f) => !f.startsWith('dist/') && (f.endsWith('.js') || f.endsWith('.ts') || f.endsWith('.sql') || f.endsWith('.py')));
  const packageFiles = ['src/package.ts', 'src/detectors.ts', 'src/package.js', 'src/detectors.js',
    'src/package-test.ts', 'src/package-test.js', 'package.js', 'detectors.js'];
  // const funcFiles = jsTsFiles.filter((f) => packageFiles.includes(f)); 
  const errors: string[] = [];
  const warnings: string[] = [];
  const packageFilePath = path.join(packagePath, 'package.json');
  const json: PackageFile = JSON.parse(fs.readFileSync(packageFilePath, { encoding: 'utf-8' }));

  const webpackConfigPath = path.join(packagePath, 'webpack.config.js');
  const isWebpack = fs.existsSync(webpackConfigPath);
  let isReleaseCandidateVersion: boolean = false;
  let externals: { [key: string]: string } | null = null;

  if (/\d+.\d+.\d+-rc(.[A-Za-z0-9]*.[A-Za-z0-9]*)?/.test(json.version))
    isReleaseCandidateVersion = true;
  if (isWebpack) {
    const content = fs.readFileSync(webpackConfigPath, { encoding: 'utf-8' });
    externals = extractExternals(content);
    if (externals)
      errors.push(...checkImportStatements(packagePath, jsTsFiles, externals));
  }
  errors.push(...checkSourceMap(packagePath));
  errors.push(...checkNpmIgnore(packagePath));
  warnings.push(...checkScriptNames(packagePath));
  const [signatureWarnings, signatureErrors] = checkFuncSignatures(packagePath, jsTsFiles);
  warnings.push(...signatureWarnings);
  errors.push(...signatureErrors);
  errors.push(...checkPackageFile(packagePath, json, { isWebpack, externals, isReleaseCandidateVersion }));
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
    return { value, message };
  }

  for (const file of files) {
    const content = fs.readFileSync(path.join(packagePath, file), { encoding: 'utf-8' });
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

export function checkFuncSignatures(packagePath: string, files: string[]): [string[], string[]] {
  const warnings: string[] = [];
  const errors: string[] = [];
  const checkFunctions: { [role: string]: FuncValidator } = {
    app: ({ name }: { name?: string }) => {
      let value = true;
      let message = '';

      if (name && typeof name === 'string') {
        const lowerCaseName = name.toLocaleLowerCase();
        if (lowerCaseName.startsWith('app')) {
          value = false;
          message += 'Prefix "App" is not needed. Consider removing it.\n';
        }
        if (lowerCaseName.endsWith('app')) {
          value = false;
          message += 'Postfix "App" is not needed. Consider removing it.\n';
        }
      }

      return { value, message };
    },
    semTypeDetector: ({ inputs, outputs }: { inputs: FuncParam[], outputs: FuncParam[] }) => {
      let value = true;
      let message = '';

      if (inputs.length !== 1 || inputs[0].type !== 'column') {
        value = false;
        message += 'Semantic type detectors must have one input of type "column"\n';
      }

      if (outputs.length !== 1 || outputs[0].type !== 'string') {
        value = false;
        message += 'Semantic type detectors must have one output of type "string"\n';
      }

      return { value, message };
    },
    cellRenderer: ({ inputs, outputs }: { inputs: FuncParam[], outputs: FuncParam[] }) => {
      let value = true;
      let message = '';

      if (inputs.length !== 0) {
        value = false;
        message += 'Cell renderer functions should take no arguments\n';
      }

      if (outputs.length !== 1 || outputs[0].type !== 'grid_cell_renderer') {
        value = false;
        message += 'Cell renderer functions must have one output of type "grid_cell_renderer"\n';
      }

      return { value, message };
    },
    viewer: ({ inputs, outputs }: { inputs: FuncParam[], outputs: FuncParam[] }) => {
      let value = true;
      let message = '';

      if (inputs.length !== 0) {
        value = false;
        message += 'Viewer functions should take no arguments\n';
      }

      if (outputs.length > 1 || (outputs.length === 1 && outputs[0].type !== 'viewer')) {
        value = false;
        message += 'Viewers must have one output of type "viewer"\n';
      }

      return { value, message };
    },
    fileViewer: ({ inputs, outputs, tags }: { inputs: FuncParam[], outputs: FuncParam[], tags?: string[] }) => {
      let value = true;
      let message = '';

      if (tags == null || (tags.length !== 1 && tags[0] !== 'fileViewer')) {
        value = false;
        message += 'File viewers must have only one tag: "fileViewer"\n';
      }

      if (inputs.length !== 1 || inputs[0].type !== 'file') {
        value = false;
        message += 'File viewers must have one input of type "file"\n';
      }

      if (outputs.length !== 1 || outputs[0].type !== 'view') {
        value = false;
        message += 'File viewers must have one output of type "view"\n';
      }

      return { value, message };
    },
    fileExporter: ({ description }: { description?: string }) => {
      let value = true;
      let message = '';

      if (description == null || description === '') {
        value = false;
        message += 'File exporters should have a description parameter\n';
      }

      return { value, message };
    },
    packageSettingsEditor: ({ outputs }: { outputs: FuncParam[] }) => {
      let value = true;
      let message = '';

      if (!(outputs.length === 1 && outputs[0].type === 'widget')) {
        value = false;
        message += 'Package settings editors must have one output of type "widget"\n';
      }

      return { value, message };
    },
    params: ({ inputs, outputs }: { inputs: FuncParam[], outputs: FuncParam[] }) => {
      let value = true;
      let message = '';

      for (const input of inputs) {
        if (!(input.name && input.type)) {
          value = false;
          message += `Function has no name or type of input parameter\n`;
        }
      }
      for (const output of outputs) {
        if (!(output.name && output.type)) {
          value = false;
          message += `Function has no name or type of output parameter\n`;
        }
      }

      return { value, message };
    },
  };
  const functionRoles = Object.keys(checkFunctions);

  for (const file of files) {
    const content = fs.readFileSync(path.join(packagePath, file), { encoding: 'utf-8' });
    const functions = getFuncMetadata(content, file.split('.').pop() ?? 'ts');

    for (const f of functions) {
      const paramsCheck = checkFunctions.params(f);
      if (!paramsCheck.value)
        warnings.push(`File ${file}, function ${f.name}:\n${paramsCheck.message}`);

      const roles = functionRoles.filter((role) => f.tags?.includes(role));
      if (roles.length > 1)
        warnings.push(`File ${file}, function ${f.name}: several function roles are used (${roles.join(', ')})`);
      else if (roles.length === 1) {
        const vr = checkFunctions[roles[0]](f);
        if (!vr.value) {
          warnings.push(`File ${file}, function ${f.name}:\n${vr.message}`);
        }
      }
      let wrongInputNames = f.inputs.filter((e) => forbiddenNames.includes(e?.name ?? ''))
      if (f.name && f.name !== 'postprocess') {
        if (namesInFiles.has(f.name))
          namesInFiles.get(f.name)?.push(file);
        else
          namesInFiles.set(f.name, [file]);
      }

      if (wrongInputNames.length > 0)
        errors.push(`File ${file}, function ${f.name}: Wrong input names: (${wrongInputNames.map((e) => e.name).join(', ')})`);
      if (f.isInvalidateOnWithoutCache)
        errors.push(`File ${file}, function ${f.name}: Can't use invalidateOn without cache, please follow this example: 'meta.cache.invalidateOn'`);

      if (f.cache)
        if (!utils.cacheValues.includes(f.cache))
          errors.push(`File ${file}, function ${f.name}: unsupposed variable for cache : ${f.cache}`);
      if (f.invalidateOn)
        if (!utils.isValidCron(f.invalidateOn))
          errors.push(`File ${file}, function ${f.name}: unsupposed variable for invalidateOn : ${f.invalidateOn}`);
    }
  }
  for(const [name, files] of namesInFiles){
    if(files.length > 1)
      errors.push(`Duplicate names ('${name}'): \n  ${files.join('\n  ')}`);
  }

  return [warnings, errors];
}

const sharedLibExternals: { [lib: string]: {} } = {
  'common/html2canvas.min.js': { 'exceljs': 'ExcelJS' },
  'common/exceljs.min.js': { 'html2canvas': 'html2canvas' },
  'common/ngl_viewer/ngl.js': { 'NGL': 'NGL' },
  'common/openchemlib-full.js': { 'openchemlib/full': 'OCL' },
  'common/codemirror/codemirror.js': { 'codemirror': 'CodeMirror' },
  'common/vue.js': { 'vue': 'Vue' }
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

    for (let dependency of Object.keys(json.dependencies ?? {})) {
      if (/\d+.\d+.\d+-rc(.[A-Za-z0-9]*.[A-Za-z0-9]*)?/.test((json.dependencies ?? {})[dependency])) {
        hasRCDependency = true;
        break;
      }
    }

    if (!hasRCDependency) {
      for (let dependency of Object.keys(json.dependencies ?? {})) {
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
    clf = fs.readFileSync(path.join(packagePath, 'CHANGELOG.md'), { encoding: 'utf-8' });
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
    const configJson: string = fs.readFileSync(tsconfigFilePath, { encoding: 'utf-8' }); // cant convert to json because file contains comments

    if (!(new RegExp('"sourceMap"\\s*:\\s*true')).test(configJson))
      warnings.push('ts config doesnt contain source map');

    const webpackConfigJson: string = fs.readFileSync(webpackConfigFilePath, { encoding: 'utf-8' }); // cant convert to json because file contains comments

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
    const npmIgnoreContent: string = fs.readFileSync(path.join(...[packagePath, '.npmignore']), { encoding: 'utf-8' });
    for (const row of npmIgnoreContent.split('\n')) {
      if ((row.match(new RegExp('\\s*dist\\/?\\s*$'))?.length || -1) > 0) {
        warnings.push('there is dist directory in .npmignore')
        break;
      }
    }
  }
  else
    warnings.push('.npmignore doesnt exists')

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
  entries.forEach(entry => {
    const entryPath = path.join(directoryPath, entry);
    const stat = fs.statSync(entryPath);

    if (stat.isFile()) {
      fileNames.push(entry);
    } else if (stat.isDirectory() && !excludedFilesToCheck.includes(entry)) {
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

function getFuncMetadata(script: string, fileExtention: string): FuncMetadata[] {
  const funcData: FuncMetadata[] = [];
  let isHeader = false;
  let data: FuncMetadata = { name: '', inputs: [], outputs: [] };

  for (const line of script.split('\n')) {
    if (!line)
      continue;
    //@ts-ignore
    const match = line.match(utils.fileParamRegex[fileExtention]);
    if (match) {
      if (!isHeader)
        isHeader = true;
      const param = match[1];
      if (param === 'name')
        data.name = line.match(utils.nameAnnRegex)?.[2]?.toLocaleLowerCase();
      else if (param === 'description')
        data.description = match[2];
      else if (param === 'input') {
        data.inputs.push({ type: match[2], name: match[3] });
      }
      else if (param === 'output')
        data.outputs.push({ type: match[2], name: match[3] });
      else if (param === 'tags') {
        data.tags = match.input && match[3] ? match.input.split(':')[1].split(',').map((t) => t.trim()) : [match[2]];
      }
      else if (param === 'meta.cache') {
        data.cache = line.split(':').pop()?.trim();
      }
      else if (param === 'meta.cache.invalidateOn') {
        data.invalidateOn = line.split(':').pop()?.trim();
      }
      else if (param === 'meta.invalidateOn') {
        data.isInvalidateOnWithoutCache = true;
      }
    }
    if (isHeader) {
      const nm = line.match(utils.nameRegex);
      if (nm)
        data.name = nm[1]?.toLocaleLowerCase();
      if (data.name && !match) {
        funcData.push(data);
        data = { name: '', inputs: [], outputs: [] };
        isHeader = false;
      }
    }
  }

  return funcData;
}

interface CheckArgs {
  _: string[],
  r?: boolean,
  recursive?: boolean,
  soft?: boolean,
}
