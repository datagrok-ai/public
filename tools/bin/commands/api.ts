import fs from 'fs';
import path from 'path';
import walk from 'ignore-walk';
import * as utils from '../utils/utils';
import * as color from '../utils/color-utils';


function generateQueryWrappers(): void {
  const curDir = process.cwd();
  const queriesDir = path.join(curDir, 'queries');
  if (!fs.existsSync(queriesDir)) {
    color.warn(`Directory ${queriesDir} not found`);
    console.log('Skipping API generation for queries...');
    return;
  }

  const packagePath = path.join(curDir, 'package.json');
  const _package = JSON.parse(fs.readFileSync(packagePath, { encoding: 'utf-8' }));

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

    const queries = script.split('--end').map((q) => q.trim()).filter((q) => q.length > 0);
    for (const q of queries) {
      const name = utils.getScriptName(q, utils.commentMap[utils.queryExtension]);
      if (!name) continue;
      let tb = new utils.TemplateBuilder(utils.queryWrapperTemplate)
        .replace('FUNC_NAME', name)
        .replace('FUNC_NAME_LOWERCASE', name)
        .replace('PACKAGE_NAMESPACE', _package.name);

      const inputs = utils.getScriptInputs(q, utils.commentMap[utils.queryExtension]);
      const outputType = utils.getScriptOutputType(q, utils.commentMap[utils.queryExtension]);
      tb.replace('PARAMS_OBJECT', inputs)
        .replace('TYPED_PARAMS', inputs)
        // The query output, if omitted, is a dataframe
        .replace('OUTPUT_TYPE', outputType === 'void' ? utils.dgToTsTypeMap['dataframe'] : outputType);
      wrappers.push(tb.build());
    }
  }

  const srcDir = path.join(curDir, 'src');
  const queryFileName = 'queries-api.ts';
  const queryFilePath = path.join(fs.existsSync(srcDir) ? srcDir : curDir, queryFileName);
  if (fs.existsSync(queryFilePath)) {
    color.warn(`The file ${queryFilePath} already exists`);
    console.log('Rewriting its contents...');
  }

  const sep = '\n';
  fs.writeFileSync(queryFilePath, utils.dgImports + sep + wrappers.join(sep.repeat(2)) + sep, 'utf8');
  color.success(`Successfully generated file ${queryFileName}${sep}`);
}

function generateScriptWrappers(): void {
  const curDir = process.cwd();
  const scriptsDir = path.join(curDir, 'scripts');
  if (!fs.existsSync(scriptsDir)) {
    color.warn(`Directory ${scriptsDir} not found`);
    console.log('Skipping API generation for scripts...');
    return;
  }

  const packagePath = path.join(curDir, 'package.json');
  const _package = JSON.parse(fs.readFileSync(packagePath, { encoding: 'utf-8' }));

  const files = walk.sync({
    path: './scripts',
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

    let tb = new utils.TemplateBuilder(utils.scriptWrapperTemplate)
      .replace('FUNC_NAME', name)
      .replace('FUNC_NAME_LOWERCASE', name)
      .replace('PACKAGE_NAMESPACE', _package.name);

    const inputs = utils.getScriptInputs(script);
    const outputType = utils.getScriptOutputType(script);
    tb.replace('PARAMS_OBJECT', inputs)
      .replace('TYPED_PARAMS', inputs)
      .replace('OUTPUT_TYPE', outputType);

    wrappers.push(tb.build(1));
  }

  const srcDir = path.join(curDir, 'src');
  const funcFileName = 'scripts-api.ts';
  const funcFilePath = path.join(fs.existsSync(srcDir) ? srcDir : curDir, funcFileName);
  if (fs.existsSync(funcFilePath)) {
    color.warn(`The file ${funcFilePath} already exists`);
    console.log('Rewriting its contents...');
  }

  const sep = '\n';
  const scriptApi = new utils.TemplateBuilder(utils.namespaceTemplate)
    .replace('PACKAGE_NAMESPACE', _package.name)
    .replace('NAME', wrappers.join(sep.repeat(2)));
  fs.writeFileSync(funcFilePath, utils.dgImports + sep + scriptApi.build() + sep, 'utf8');
  color.success(`Successfully generated file ${funcFileName}${sep}`);
}

export function api(args: { _: string[] }): boolean {
  const nOptions = Object.keys(args).length - 1;
  if (args['_'].length !== 1 || nOptions > 0) return false;
  if (!utils.isPackageDir(process.cwd())) {
    color.error('File `package.json` not found. Run the command from the package directory');
    return false;
  }
  generateScriptWrappers();
  generateQueryWrappers();
  return true;
}
