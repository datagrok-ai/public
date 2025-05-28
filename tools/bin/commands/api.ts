import fs from 'fs';
import path from 'path';
import walk from 'ignore-walk';
import * as utils from '../utils/utils';
import * as color from '../utils/color-utils';

const sep = '\n';

const packageFuncDirs = ['package.ts', 'package.g.ts'];
const apiFile = 'package-api.ts';

const curDir = process.cwd();
const srcDir = path.join(curDir, 'src');
const funcFilePath = path.join(fs.existsSync(srcDir) ? srcDir : curDir, apiFile);
const packagePath = path.join(curDir, 'package.json');
const names = new Set<string>();

const _package = JSON.parse(fs.readFileSync(packagePath, { encoding: 'utf-8' }));

function generateQueryWrappers(): void {
  const queriesDir = path.join(curDir, 'queries');
  if (!fs.existsSync(queriesDir)) {
    color.warn(`Directory ${queriesDir} not found`);
    console.log('Skipping API generation for queries...');
    return;
  }

  const files = walk.sync({
    path: './queries',
    ignoreFiles: ['.npmignore', '.gitignore'],
  });

  const wrappers = [];
  for (const file of files) {
    if (!file.endsWith(utils.queryExtension)) continue;
    const filepath = path.join(queriesDir, file);
    const script = fs.readFileSync(filepath, 'utf8');
    if (!script) continue;

    const queries = script.split(/--\s*end/).map((q) => q.trim()).filter((q) => q.length > 0);
    for (const q of queries) {
      const name = utils.getScriptName(q, utils.commentMap[utils.queryExtension]);
      if (!name) continue;

      checkNameColision(name);

      const tb = new utils.TemplateBuilder(utils.queryWrapperTemplate)
        .replace('FUNC_NAME', name)
        .replace('FUNC_NAME_LOWERCASE', name)
        .replace('PACKAGE_NAMESPACE', _package.name);

      const description = utils.getScriptDescription(q, utils.commentMap[utils.queryExtension]);
      const inputs = utils.getScriptInputs(q, utils.commentMap[utils.queryExtension]);
      const outputType = utils.getScriptOutputType(q, utils.commentMap[utils.queryExtension]);
      tb.replace('PARAMS_OBJECT', inputs)
        .replace('TYPED_PARAMS', inputs)
        .replace('FUNC_DESCRIPTION', description)
        .replace('OUTPUT_TYPE', outputType === 'void' ? utils.dgToTsTypeMap['dataframe'] : outputType);
      wrappers.push(tb.build(1));
    }
  }

  saveWrappersToFile('queries', wrappers);
}

function generateScriptWrappers(): void {
  const scriptsDir = path.join(curDir, 'scripts');
  if (!fs.existsSync(scriptsDir)) {
    color.warn(`Directory ${scriptsDir} not found`);
    console.log('Skipping API generation for scripts...');
    return;
  }

  const files = walk.sync({
    path: scriptsDir,
    ignoreFiles: ['.npmignore', '.gitignore'],
  });
  const wrappers = [];
  for (const file of files) {
    let extension: string;
    if (!utils.scriptExtensions.some((ext) => (extension = ext, file.endsWith(ext)))) continue;

    const filepath = path.join(scriptsDir, file);
    const script = fs.readFileSync(filepath, 'utf8');
    if (!script) continue;

    const name = utils.getScriptName(script, utils.commentMap[extension!]);
    if (!name) continue;
    const description = utils.getScriptDescription(script);

    checkNameColision(name)

    const tb = new utils.TemplateBuilder(utils.scriptWrapperTemplate)
      .replace('FUNC_NAME', name)
      .replace('FUNC_NAME_LOWERCASE', name)
      .replace('PACKAGE_NAMESPACE', _package.name);

    const inputs = utils.getScriptInputs(script);
    const outputType = utils.getScriptOutputType(script);
    tb.replace('PARAMS_OBJECT', inputs)
      .replace('TYPED_PARAMS', inputs)
      .replace('FUNC_DESCRIPTION', description)
      .replace('OUTPUT_TYPE', outputType);
    wrappers.push(tb.build(1));
  }

  saveWrappersToFile('scripts', wrappers);
}

function generateFunctionWrappers(): void {
  let filesToParse = packageFuncDirs.map((e) => path.join(curDir, 'src', e)).filter((e) => fs.existsSync(e));

  const annotaionRegex = /(?:\/\/[^\n]*\n)+export[^{]*/g;
  const nameRegex = /\s*export(?:\sasync)?\s*function\s*([^\s(]*)/;

  const wrappers = [];
  for (const file of filesToParse) {
    const fileData = fs.readFileSync(file, 'utf8');
    const annotations = fileData.matchAll(annotaionRegex);
    if (annotations === null)
      return;

    for (let annotation of annotations) {
      const name = (annotation[0].match(nameRegex) ?? [undefined, undefined])[1];
      const description = utils.getScriptDescription(annotation[0], utils.commentMap[utils.jsExtention]);
      if (!name)
        continue;

      const annotationInputs = utils.getScriptInputs(annotation[0], utils.commentMap[utils.jsExtention]);
      const annotationOutputDir = utils.getScriptOutputType(annotation[0], utils.commentMap[utils.jsExtention]);
      let outputType = '';

      for (let outputAnnotation of annotationOutputDir) {
        if (outputType != '') {
          outputType = 'any';
          break;
        }
        outputType = utils.dgToTsTypeMap[outputAnnotation[1]] ?? 'any';
      }

      checkNameColision(name);

      const tb = new utils.TemplateBuilder(utils.scriptWrapperTemplate)
        .replace('FUNC_NAME', name)
        .replace('FUNC_NAME_LOWERCASE', name)
        .replace('PACKAGE_NAMESPACE', _package.name)
        .replace('PARAMS_OBJECT', annotationInputs)
        .replace('FUNC_DESCRIPTION', description)
        .replace('TYPED_PARAMS', annotationInputs)
        .replace('OUTPUT_TYPE', outputType);
      wrappers.push(tb.build(1));
    }
  }

  saveWrappersToFile('funcs', wrappers);
}

function saveWrappersToFile(namespaceName: string, wrappers: string[]) {
  if (!fs.existsSync(funcFilePath))
    createApiFile()

  const scriptApi = new utils.TemplateBuilder(utils.namespaceTemplate)
    .replace('PACKAGE_NAMESPACE', namespaceName)
    .replace('NAME', wrappers.join(sep.repeat(2)));
  fs.appendFileSync(funcFilePath, sep + scriptApi.build() + sep);
  color.success(`Successfully generated file ${apiFile}${sep}`);
}

function createApiFile() {
  if (fs.existsSync(funcFilePath)) {
    color.warn(`The file ${funcFilePath} already exists`);
    console.log('Rewriting its contents...');
  }
  fs.writeFileSync(funcFilePath, utils.dgImports + sep, 'utf8');
}

function checkNameColision(name: string) {
  if (names.has(name))
    console.log('There is collision in name ' + name)
  names.add(name);
}

export function api(args: { _: string[] }): boolean {
  const nOptions = Object.keys(args).length - 1;
  if (args['_'].length !== 1 || nOptions > 0) return false;
  if (!utils.isPackageDir(process.cwd())) {
    color.error('File `package.json` not found. Run the command from the package directory');
    return false;
  }
  createApiFile();
  generateScriptWrappers();
  generateQueryWrappers();
  generateFunctionWrappers();
  return true;
}
