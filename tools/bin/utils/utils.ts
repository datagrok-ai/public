import fs from 'fs';
import path, {sep} from 'path';
import {exec} from 'child_process';
import {promisify} from 'util';

const execAsync = promisify(exec);


/* Waits [ms] milliseconds */
export async function delay(ms: number) {
  await new Promise((r) => setTimeout(r, ms));
}

export function isEmpty(dir: string): boolean {
  return fs.readdirSync(dir).length === 0;
}

export function isPackageDir(dir: string): boolean {
  return fs.existsSync(path.join(dir, 'package.json'));
}

export function kebabToCamelCase(s: string, firstUpper: boolean = true): string {
  s = s.replace(/-./g, (x) => x.toUpperCase()[1]);
  return (firstUpper ? s[0].toUpperCase() : s[0].toLowerCase()) + s.slice(1);
}

export function descriptionToComment(s: string) {
  if (s.length === 0)
    return '';
  return '/**\n' + s + '\n*/\n';
}

export function spaceToCamelCase(s: string, firstUpper: boolean = true): string {
  s = s.replace(/\s+./g, (x) => x[x.length - 1].toUpperCase());
  return (firstUpper ? s[0].toUpperCase() : s[0].toLowerCase()) + s.slice(1);
}

export function wordsToCamelCase(s: string, firstUpper: boolean = true): string {
  const m = s.match(/^[A-Z]+/); // only abbreviation handling
  const lastIndex = m ? m[0].length : 1;
  return s.slice(0, lastIndex).toLowerCase() + s.slice(lastIndex);
}

export function camelCaseToKebab(s: string): string {
  return s.replace(/[A-Z]/g, (char: string, index: number) => index == 0 ? char.toLowerCase() : '-' + char.toLowerCase());
}

export function removeScope(name: string): string {
  const split = name.split('/');
  return split[split.length - 1];
}

export function mapURL(conf: Config): Indexable {
  const urls: Indexable = {};
  for (const server in conf['servers'])
    urls[conf['servers'][server]['url']] = server;

  return urls;
}

export function friendlyNameToName(s: string, firstUpper: boolean = true): string {
  let out = '';
  let cap = true;
  let firstWordUpperCase = false;
  let start = true;
  s = (s ?? '').trim();
  const letterRegex = /[A-Za-z]/;
  const digitRegex = /\d/;
  const isUpper = (s: string) => /[A-Z]/.test(s);
  for (let i = 0; i < s.length; i++) {
    if (!letterRegex.test(s[i]) && !digitRegex.test(s[i]))
      cap = true;
    else {
      if (start && digitRegex.test(s[i]))
        continue;
      firstWordUpperCase = start ? isUpper(s[i]) : firstWordUpperCase && isUpper(s[i]);
      out += !firstUpper && (start || firstWordUpperCase && !cap) ? s[i].toLowerCase() : cap ? s[i].toUpperCase() : s[i];
      cap = false;
      start = false;
    }
  }
  return out;
}

export const replacers: Indexable = {
  NAME: (s: string, name: any) => s.replace(/#{NAME}/g, name),
  NAME_TITLECASE: (s: string, name: string) => s.replace(/#{NAME_TITLECASE}/g, name[0].toUpperCase() + name.slice(1).toLowerCase()),
  NAME_LOWERCASE: (s: string, name: string) => s.replace(/#{NAME_LOWERCASE}/g, name.toLowerCase()),
  NAME_PREFIX: (s: string, name: string) => s.replace(/#{NAME_PREFIX}/g, name.slice(0, 3)),
  PACKAGE_DETECTORS_NAME: (s: string, name: string) => s.replace(/#{PACKAGE_DETECTORS_NAME}/g, kebabToCamelCase(name)),
  PACKAGE_NAMESPACE: (s: string, name: string) => s.replace(/#{PACKAGE_NAMESPACE}/g, name),
  FUNC_DESCRIPTION: (s: string, desc: string) => s.replace(/#{FUNC_DESCRIPTION}/g, descriptionToComment(desc)),
  FUNC_NAME: (s: string, name: string) => s.replace(/#{FUNC_NAME}/g, friendlyNameToName(name)),
  FUNC_NAME_LOWERCASE: (s: string, name: string) => s.replace(/#{FUNC_NAME_LOWERCASE}/g, friendlyNameToName(name, false)),
  PARAMS_OBJECT: (s: string, params: { name?: string; type?: string }[]) => s.replace(/#{PARAMS_OBJECT}/g, params.length ?
    `{ ${params.map((p) => p.name).join(', ')} }` : `{}`),
  OUTPUT_TYPE: (s: string, type: string) => s.replace(/#{OUTPUT_TYPE}/g, type),
  TYPED_PARAMS: (s: string, params: { name?: string; type?: string; isOptional?: boolean; undefinable: boolean, nullable?: boolean }[]) => s.replace(/#{TYPED_PARAMS}/g,
    params.map((p) => `${p.name}${p.isOptional ? '?' : ''}: ${p.type} ${p.undefinable ? '| undefined' : '' }${p.nullable ? '| null' : ''}`).join(', ')),
};


export class TemplateBuilder {
  static sep = '\n';
  static indentSize = 2;
  template: string;
  constructor(template: string) {
    this.template = template;
  }

  replace(pattern: string, value: string | object[]) {
    this.template = replacers[pattern](this.template, value);
    return this;
  }

  build(indent: number = 0) {
    if (indent) {
      this.template = this.template
        .split(TemplateBuilder.sep)
        .map((line) => ' '.repeat(indent * TemplateBuilder.indentSize) + line)
        .join(TemplateBuilder.sep);
    }
    return this.template;
  }
}

export const scriptLangExtMap: Indexable = {
  javascript: 'js',
  julia: 'jl',
  node: 'js',
  octave: 'm',
  python: 'py',
  r: 'R',
};

export const commentMap: Indexable = {
  '.js': '//',
  '.jl': '#',
  '.m': '#',
  '.py': '#',
  '.R': '#',
  '.sql': '--',
};

export const queryExtension = '.sql';
export const jsExtention = '.js';
export const scriptExtensions = ['.jl', '.m', '.py', '.R'];
export function checkScriptLocation(filepath: string): boolean {
  if (!(filepath.startsWith('scripts/') || filepath.startsWith('projects/') || filepath.startsWith('dockerfiles/') || filepath.startsWith('python/')) &&
    scriptExtensions.some((ext: any) => filepath.endsWith(ext)))
    return false;

  return true;
};

export function getScriptName(script: string, comment: string = '#'): string | null {
  const regex = new RegExp(`${comment}\\s*name:\\s*(.*)`);
  const match = script.match(regex);
  return match ? match[1]?.trim() : null;
};

export function getParam(name: string, script: string, comment: string = '#'): string | null {
  const regex = new RegExp(`${comment}\\s*${name}:\\s*(.*)`);
  const match = script.match(regex);
  return match ? match[1]?.trim() : null;
};

export const cacheValues = ['all', 'server', 'client', 'true'];

export function isValidCron(cronExpression: string): boolean {
  const cronRegex = /^(\*|([0-9]|1[0-9]|2[0-9]|3[0-9]|4[0-9]|5[0-9])|\*\/([0-9]|1[0-9]|2[0-9]|3[0-9]|4[0-9]|5[0-9])) (\*|([0-9]|1[0-9]|2[0-3])|\*\/([0-9]|1[0-9]|2[0-3])) (\*|([1-9]|1[0-9]|2[0-9]|3[0-1])|\*\/([1-9]|1[0-9]|2[0-9]|3[0-1])) (\*|([1-9]|1[0-2])|\*\/([1-9]|1[0-2])) (\*|([0-6])|\*\/([0-6]))$/;
  return cronRegex.test(cronExpression);
}

export const dgToTsTypeMap: Indexable = {
  int: 'number',
  double: 'number',
  bigint: 'bigint',
  bool: 'boolean',
  string: 'string',
  dataframe: 'DG.DataFrame',
  column: 'DG.Column',
  column_list: 'string[]',
  file: 'DG.FileInfo',
  view: 'DG.View',
  void: 'void',
};

export const propertyTypes = [
  'bool', 'int', 'double', 'string', 'datetime', 'object',
  'column', 'dataframe', 'bitset', 'cell', 'string_list', 'map',
];

export const headerTags = [
  'name', 'description', 'help-url', 'input', 'output', 'tags',
  'sample', 'language', 'returns', 'test', 'sidebar', 'condition',
  'top-menu', 'environment', 'require', 'editor-for', 'schedule',
  'reference', 'editor', 'meta', 'connection', 'friendlyName',
];

export const fileParamRegex = {
  py: new RegExp(`^#\\s*(?!\\s*#)([^:]+):\\s+([^\\s\\[\\{]+) ?([^\\s\\[\\{]+)?`),
  ts: new RegExp(`^//\\s*(?!\\s*//)([^:]+):\\s+([^\\s\\[\\{]+) ?([^\\s\\[\\{]+)?`),
  js: new RegExp(`^//\\s*(?!\\s*//)([^:]+):\\s+([^\\s\\[\\{]+) ?([^\\s\\[\\{]+)?`),
  sql: new RegExp(`^--\\s*([^:]+):\\s+([^\\s\\[\\{]+) ?([^\\s\\[\\{]+)?`),
};

export const nameAnnRegex = /\s*(name[^:]*): ([^\n\r\[\{]+)/;

export const nameRegex = /(?:|(?:static)(?:export )(?:async ))\s+function\s+([a-zA-Z_][a-zA-Z0-9_$]*)\s*\((.*?).*/;

export const absUrlRegex = new RegExp('^(?:[a-z+]+:)?//', 'i');

export function getScriptOutputType(script: string, comment: string = '#'): string {
  const regex = /\s*output:\s?([a-z_]+)\s*([^\s]*)/g;
  const matches = script.matchAll(regex);
  if (!matches) return 'void';
  let resType = 'void';
  let firstItemName = '';
  let wasSecond = false;
  for (const match of matches) {
    if (resType === 'void') {
      resType = dgToTsTypeMap[match[1]] ?? 'any';
      firstItemName = match[2];
    } else {
      if (!wasSecond) {
        resType = `${firstItemName}: ${resType}`;
        wasSecond = true;
      }
      resType = [resType, `${match[2]}: ${dgToTsTypeMap[match[1]] ?? 'any'}`].join(', ');
    }
  }
  return wasSecond ? `{${resType}}` : resType;
};

export function getScriptInputs(script: string, comment: string = '#'): object[] {
  const regex = new RegExp(`${comment}\\s*input:\\s?([a-z_]+)(?:<[^>]*>)?\\s+(\\w+)(?:[^{\\n]*{[^}\\n]*})?`, 'g');
  const testOptional = (inputAnnotation: string) => /isOptional\s*:\s*true/.test(inputAnnotation) || /optional\s*:\s*true/.test(inputAnnotation);
  const inputAnnotations = [...script.matchAll(regex)];
  let firstTsValidOptionalIdx: number = inputAnnotations.length-1;
  for (; firstTsValidOptionalIdx >= 0; firstTsValidOptionalIdx--) {
    const annotation = inputAnnotations[firstTsValidOptionalIdx][0];
    if (!testOptional(annotation))
      break;
  }
  const inputs = [];
  for (const [idx, match] of inputAnnotations.entries()) {
    const hasOptionalAnnotation = testOptional(match[0]);
    const isOptional = hasOptionalAnnotation && idx >= firstTsValidOptionalIdx;
    const undefinable = hasOptionalAnnotation && idx < firstTsValidOptionalIdx;
    const nullable = /nullable\s*:\s*true/.test(match[0]);
    const type = dgToTsTypeMap[match[1]] || 'any';
    const name = match[2];
    inputs.push({type, name, isOptional, undefinable, nullable});
  }
  return inputs;
};

export function getScriptDescription(script: string, comment: string = '#'): string {
  const regex = new RegExp(`${comment}\\s*description:\\s([^\n]*)`);
  const rexegRes = script.match(regex) || [];
  const desc = rexegRes[1] || '';
  return desc;
};

export const dgImports = `import * as grok from 'datagrok-api/grok';\nimport * as DG from 'datagrok-api/dg';\n`;

export const scriptWrapperTemplate = `#{FUNC_DESCRIPTION}export async function #{FUNC_NAME_LOWERCASE}(#{TYPED_PARAMS}): Promise<#{OUTPUT_TYPE}> {
  return await grok.functions.call('#{PACKAGE_NAMESPACE}:#{FUNC_NAME}', #{PARAMS_OBJECT});
}`;

export const queryWrapperTemplate = `#{FUNC_DESCRIPTION}export async function #{FUNC_NAME_LOWERCASE}(#{TYPED_PARAMS}): Promise<#{OUTPUT_TYPE}> {
  return await grok.data.query('#{PACKAGE_NAMESPACE}:#{FUNC_NAME}', #{PARAMS_OBJECT});
}`;

export const namespaceTemplate = `export namespace #{PACKAGE_NAMESPACE} {\n#{NAME}\n}`;

export interface Config {
  servers: {
    [alias: string]: {
      url: string,
      key: string
    }
  },
  default: string,
}

export interface Indexable { [key: string]: any }

export async function runScript(script: string, path: string, verbose: boolean = false) {
  try {
    const {stdout, stderr} = await execAsync(script, {cwd: path});
    if (stderr && verbose) 
      console.error(`Warning/Error: ${stderr}`);
    
    if (stdout && verbose) 
      console.log(`Output: ${stdout}`);
    
  } catch (error: any) {
    console.error(`Execution failed: ${error.message}`);
    throw new Error(`Error executing '${script}'. Error message: ${error.message}`);
  }
}

export function setHost(host: any, configFile: any) {
  if (host) {
    if (host in configFile.servers) {
      process.env.HOST = host;
      console.log('Environment variable `HOST` is set to', host);
    } else {
      console.error(`Unknown server alias. Please add it to Config File`);
      return false;
    }
  } else if (configFile.default) {
    process.env.HOST = configFile.default;
    console.log('Environment variable `HOST` is set to', configFile.default);
  }
}

export async function runAll(packagesDir: string, command: string, options: Record<string, any>, packagesToLoad:string[] = ['all']) {
  const packages = fs.readdirSync(packagesDir);
  const commandToRun = `${command} ${optionsToString(options)}`;
  for (const packageName of packages) {
    const packagePath = path.join(...[packagesDir, packageName]);
    const packageJsonPath = path.join(...[packagesDir, packageName, 'package.json']);
    if (fs.statSync(packagePath).isDirectory() && fs.existsSync(packageJsonPath)) {
      try {
        console.log(`${packagePath}: ${commandToRun} - Started`);
        await runScript(commandToRun, packagePath, options.verbose);
        console.log(`${packagePath}: ${commandToRun} - Finished`);
      } catch (e) {
        console.log(`${packagePath}: ${commandToRun} - Error`);
      }
    }
  }
  return;
}

function optionsToString(options: Record<string, any>) {
  const parts: string[] = []; 

  for (const [key, value] of Object.entries(options)) {
    if (key === '_' || key === 'all') continue;

    if (typeof value === 'boolean') {
      if (value) 
        parts.push(`--${key}`);
      
    } else if (value !== undefined && value !== null) 
      parts.push(`--${key}="${value}"`);
    
  }

  return parts.join(' ');
}
