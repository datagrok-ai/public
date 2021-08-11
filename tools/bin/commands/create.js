const fs = require('fs');
const path = require('path');
const os = require('os');
const yaml = require('js-yaml');
const help = require('../ent-helpers.js');
const utils = require('../utils.js');
var exec = require('child_process').exec;

module.exports = {
  create: create
};

const platform = os.platform();
const curDir = process.cwd();
const curFolder = path.basename(curDir);

const grokDir = path.join(os.homedir(), '.grok');
const confPath = path.join(grokDir, 'config.yaml');

const templateDir = path.join(path.dirname(path.dirname(__dirname)), 'package-template');
const confTemplateDir = path.join(path.dirname(path.dirname(__dirname)), 'config-template.yaml');

const confTemplate = yaml.safeLoad(fs.readFileSync(confTemplateDir));

function createDirectoryContents(name, config, templateDir, packageDir, ide = '', ts = false) {
  const filesToCreate = fs.readdirSync(templateDir);

  filesToCreate.forEach(file => {
    const origFilePath = path.join(templateDir, file);
    let copyFilePath = path.join(packageDir, file);
    const stats = fs.statSync(origFilePath);
    if (stats.isFile()) {
      if (file === 'package.png') {
        fs.writeFileSync(copyFilePath, fs.readFileSync(origFilePath, 'base64'), 'base64');
        return false;
      }
      let contents = fs.readFileSync((file === 'webpack.config.js' && ts) ?
        path.join(templateDir, 'ts.webpack.config.js') : origFilePath, 'utf8');
      contents = contents.replace(/#{PACKAGE_NAME}/g, name);
      contents = contents.replace(/#{PACKAGE_DETECTORS_NAME}/g, utils.kebabToCamelCase(name));
      contents = contents.replace(/#{PACKAGE_NAME_LOWERCASE}/g, name.toLowerCase());
      contents = contents.replace(/#{PACKAGE_NAME_LOWERCASE_WORD}/g, name.replace(/-/g, '').toLowerCase());
      contents = contents.replace(/#{GROK_HOST_ALIAS}/g, config.default);
      contents = contents.replace(/#{GROK_HOST}/g, /localhost|127\.0\.0\.1/.test(config['servers'][config.default]['url']) ?
        'http://localhost:63343/login.html' : (new URL(config['servers'][config.default]['url'])).origin);
      if (file === 'package.json') {
        // Generate scripts for non-default servers from `config.yaml`
        let package = JSON.parse(contents);
        for (let server in config.servers) {
          if (server === config.default) continue;
          package['scripts'][`debug-${name.toLowerCase()}-${server}`] = `grok publish ${server} --rebuild`;
          package['scripts'][`release-${name.toLowerCase()}-${server}`] = `grok publish ${server} --rebuild --release`;
        }
        if (ts) Object.assign(package.dependencies, { 'ts-loader': 'latest', 'typescript': 'latest' });
        contents = JSON.stringify(package, null, '\t');
      }
      if (file === 'package.js' && ts) copyFilePath = path.join(packageDir, 'package.ts');
      if (file === 'tsconfig.json' && !ts) return false;
      if (file === 'ts.webpack.config.js') return false;
      if (file === 'gitignore') {
        copyFilePath = path.join(packageDir, '.gitignore');
        if (ts) contents += '\n# Emitted *.js files\nsrc/**/*.js\n';
      }
      fs.writeFileSync(copyFilePath, contents, 'utf8');
    } else if (stats.isDirectory()) {
      if (file === '.vscode' && !(ide == 'vscode' && platform == 'win32')) return;
      fs.mkdirSync(copyFilePath);
      // recursive call
      createDirectoryContents(name, config, origFilePath, copyFilePath, ide, ts);
    }
  })
}

function create(args) {
  const nOptions = Object.keys(args).length - 1;
  const nArgs = args['_'].length;
  if (nArgs > 2 || nOptions > 2) return false;
  if (nOptions && !Object.keys(args).slice(1).every(op => op == 'ide' || op == 'ts')) return false;

  // Create `config.yaml` if it doesn't exist yet
  if (!fs.existsSync(grokDir)) fs.mkdirSync(grokDir);
  if (!fs.existsSync(confPath)) fs.writeFileSync(confPath, yaml.safeDump(confTemplate));

  const config = yaml.safeLoad(fs.readFileSync(confPath));

  const name = nArgs === 2 ? args['_'][1] : curFolder;
  const validName = /^([A-Za-z\-_\d])+$/.test(name);

  if (validName) {
    let packageDir = curDir;
    if (curFolder !== name) {
      packageDir = path.join(packageDir, name);
      if (!fs.existsSync(packageDir)) {
        fs.mkdirSync(packageDir);
      }
    }
    if (!utils.isEmpty(packageDir)) {
      console.log();
      console.log('The package directory should be empty');
      return false;
    }
    createDirectoryContents(name, config, templateDir, packageDir, args.ide, args.ts);
    console.log(help.package(name, args.ts));
    console.log('\nRunning `npm install` to get the required dependencies...\n');
    exec('npm install', { cwd: packageDir }, (err, stdout, stderr) => {
      if (err) throw err;
      else console.log(stderr, stdout);
    });
  } else {
    console.log('Package name may only include letters, numbers, underscores, or hyphens');
  }
  return true;
}
