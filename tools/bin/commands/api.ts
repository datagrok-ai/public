import fs from 'fs';
import path from 'path';
import walk from 'ignore-walk';
import * as utils from '../utils.js';


function generateQueryWrappers(): void {
  const curDir = process.cwd();
  const queriesDir = path.join(curDir, 'queries');
  if (!fs.existsSync(queriesDir)) {
    console.log(`Directory ${queriesDir} not found\nSkipping API generation for queries...`);
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

    const name = utils.getScriptName(script, utils.commentMap[utils.queryExtension]);
    if (!name) continue;
  }
}

function generateScriptWrappers(): void {
  const curDir = process.cwd();
  const scriptsDir = path.join(curDir, 'scripts');
  if (!fs.existsSync(scriptsDir)) {
    console.log(`Directory ${scriptsDir} not found\nSkipping API generation for scripts...`);
    return;
  }

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
      .replace('FUNC_NAME_LOWERCASE', name);

    const packagePath = path.join(curDir, 'package.json');
    const _package = JSON.parse(fs.readFileSync(packagePath, { encoding: 'utf-8' }));
    tb.replace('PACKAGE_NAMESPACE', _package.name);

    const inputs = utils.getScriptInputs(script);
    const outputType = utils.getScriptOutputType(script);
    tb.replace('PARAMS_OBJECT', inputs)
      .replace('TYPED_PARAMS', inputs)
      .replace('OUTPUT_TYPE', outputType);

    wrappers.push(tb.build());
  }

  const srcDir = path.join(curDir, 'src');
  let funcFilePath = path.join(fs.existsSync(srcDir) ? srcDir : curDir, 'scripts-api.ts');
  if (fs.existsSync(funcFilePath)) {
    console.log(`The file ${funcFilePath} already exists\nRewriting its contents...`);
  }

  const sep = '\n';
  fs.writeFileSync(funcFilePath, utils.dgImports + sep + wrappers.join(sep.repeat(2)) + sep, 'utf8');
}

export function api(args: { _: string[] }): boolean {
  const nOptions = Object.keys(args).length - 1;
  if (args['_'].length !== 1 || nOptions > 0) return false;
  if (!utils.isPackageDir(process.cwd())) {
    console.log('File `package.json` not found. Run the command from the package directory');
    return false;
  }
  generateScriptWrappers();
  generateQueryWrappers();
  return true;
}
