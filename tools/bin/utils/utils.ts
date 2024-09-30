import fs from 'fs';
import path from 'path';
import { exec } from 'child_process';
import { promisify } from 'util';

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
  PACKAGE_NAMESPACE: (s: string, name: string) => s.replace(/#{PACKAGE_NAMESPACE}/g, kebabToCamelCase(removeScope(name))),
  FUNC_NAME: (s: string, name: string) => s.replace(/#{FUNC_NAME}/g, friendlyNameToName(name)),
  FUNC_NAME_LOWERCASE: (s: string, name: string) => s.replace(/#{FUNC_NAME_LOWERCASE}/g, friendlyNameToName(name, false)),
  PARAMS_OBJECT: (s: string, params: { name?: string; type?: string }[]) => s.replace(/#{PARAMS_OBJECT}/g, params.length ?
    `{ ${params.map((p) => p.name).join(', ')} }` : `{}`),
  OUTPUT_TYPE: (s: string, type: string) => s.replace(/#{OUTPUT_TYPE}/g, type),
  TYPED_PARAMS: (s: string, params: { name?: string; type?: string }[]) => s.replace(/#{TYPED_PARAMS}/g,
    params.map((p) => `${p.name}: ${p.type}`).join(', ')),
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
export const scriptExtensions = ['.jl', '.m', '.py', '.R'];
export function checkScriptLocation(filepath: string): boolean {
  if (!(filepath.startsWith('scripts/') || filepath.startsWith('projects/') || filepath.startsWith('dockerfiles/')) &&
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
};

export const propertyTypes = [
  'bool', 'int', 'double', 'string', 'datetime', 'object',
  'column', 'dataframe', 'bitset', 'cell', 'string_list', 'map',
];

export const headerTags = [
  'name', 'description', 'help-url', 'input', 'output', 'tags',
  'sample', 'language', 'returns', 'test', 'sidebar', 'condition',
  'top-menu', 'environment', 'require', 'editor-for', 'schedule',
  'reference', 'editor', 'meta',
];

export const paramRegex = new RegExp(`\/\/\\s*(${headerTags.join('|')}\\.[^:]*): *([^\\s\\[\\{]+) ?([^\\s\\[\\{]+)?`);

export const nameAnnRegex = /\/\/\s*(name[^:]*): ([^\n\r\[\{]+)/;

export const nameRegex = /(?:|static|export\s+function|export\s+async\s+function)\s+([a-zA-Z_][a-zA-Z0-9_$]*)\s*\((.*?)\).*/;

export const absUrlRegex = new RegExp('^(?:[a-z+]+:)?//', 'i');

export function getScriptOutputType(script: string, comment: string = '#'): string {
  const regex = new RegExp(`${comment}\\s*output:\\s?([a-z_]+)\\s*`);
  const match = script.match(regex);
  if (!match) return 'void';
  return dgToTsTypeMap[match[1]] || 'any';
};

export function getScriptInputs(script: string, comment: string = '#'): object[] {
  const regex = new RegExp(`${comment}\\s*input:\\s?([a-z_]+)\\s+(\\w+)`, 'g');
  const inputs = [];
  for (const match of script.matchAll(regex)) {
    const type = dgToTsTypeMap[match[1]] || 'any';
    const name = match[2];
    inputs.push({ type, name });
  }
  return inputs;
};

export const dgImports = `import * as grok from 'datagrok-api/grok';\nimport * as DG from 'datagrok-api/dg';\n\n`;

export const scriptWrapperTemplate = `export async function #{FUNC_NAME_LOWERCASE}(#{TYPED_PARAMS}): Promise<#{OUTPUT_TYPE}> {
  return await grok.functions.call('#{PACKAGE_NAMESPACE}:#{FUNC_NAME}', #{PARAMS_OBJECT});
}`;

export const queryWrapperTemplate = `export async function #{FUNC_NAME_LOWERCASE}(#{TYPED_PARAMS}): Promise<#{OUTPUT_TYPE}> {
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
    const { stdout, stderr } = await execAsync(script, { cwd: path });
    if (stderr && verbose) {
      console.error(`Warning/Error: ${stderr}`);
    }
    if (stdout && verbose) {
      console.log(`Output: ${stdout}`);
    }
  } catch (error: any) {
    console.error(`Execution failed: ${error.message}`);
    throw new Error("Cant run script");
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