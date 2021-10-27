const fs = require('fs');
const path = require('path');
const utils = require('../utils.js');
const walk = require('ignore-walk');


function generateWrappers() {
  const curDir = process.cwd();
  const scriptsDir = path.join(curDir, 'scripts');
  if (!fs.existsSync(scriptsDir)) {
    console.log(`Directory ${scriptsDir} not found\n` +
      'Place your scripts there before running this command');
    return;
  }

  const files = walk.sync({
    path: './scripts',
    ignoreFiles: ['.npmignore', '.gitignore'],
  });

  const wrappers = [];
  for (const file of files) {
    let extension = null;
    if (!utils.scriptExtensions.some((ext) => (extension = ext, file.endsWith(ext)))) continue;

    const filepath = path.join(scriptsDir, file);
    const script = fs.readFileSync(filepath, 'utf8');
    if (!script) continue;

    const name = utils.getScriptName(script);
    if (!name) continue;

    let tb = new utils.TemplateBuilder(utils.scriptWrapperTemplate)
      .replace('FUNC_NAME', name)
      .replace('FUNC_NAME_LOWERCASE', name);

    const packagePath = path.join(curDir, 'package.json');
    const package = JSON.parse(fs.readFileSync(packagePath));
    tb.replace('PACKAGE_NAMESPACE', package.name);

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

function api(args) {
  const nOptions = Object.keys(args).length - 1;
  if (args['_'].length !== 1 || nOptions > 0) return false;
  generateWrappers();
  return true;
}

module.exports = {
  api: api,
};
