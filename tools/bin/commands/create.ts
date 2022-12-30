import fs from 'fs';
import path from 'path';
import os from 'os';
import yaml from 'js-yaml';
import { exec } from 'child_process';
import { help } from '../utils/ent-helpers';
import * as utils from '../utils/utils';
import * as color from '../utils/color-utils';
import { validateConf } from '../validators/config-validator';


const platform = os.platform();
const curDir = process.cwd();
const curFolder = path.basename(curDir);

const grokDir = path.join(os.homedir(), '.grok');
const confPath = path.join(grokDir, 'config.yaml');

const templateDir = path.join(path.dirname(path.dirname(__dirname)), 'package-template');
const confTemplateDir = path.join(path.dirname(path.dirname(__dirname)), 'config-template.yaml');
const confTemplate = yaml.load(fs.readFileSync(confTemplateDir, { encoding: 'utf-8' }));

let dependencies: string[] = [];

function createDirectoryContents(name: string, config: utils.Config, templateDir: string, packageDir: string, ide: string = '', ts: boolean = true, eslint: boolean = false, jest: boolean = false) {

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
      contents = utils.replacers['PACKAGE_NAMESPACE'](contents, name);
      contents = contents.replace(/#{GROK_HOST_ALIAS}/g, config.default);
      contents = contents.replace(/#{GROK_HOST}/g, /localhost|127\.0\.0\.1/.test(config['servers'][config.default]['url']) ?
        'http://localhost:63343/login.html' : (new URL(config['servers'][config.default]['url'])).origin);
      if (file === 'package.json') {
        // Generate scripts for non-default servers from `config.yaml`
        let _package = JSON.parse(contents);
        for (let server in config.servers) {
          if (server === config.default) continue;
          _package['scripts'][`debug-${name.toLowerCase()}-${server}`] = `webpack && grok publish ${server}`;
          _package['scripts'][`release-${name.toLowerCase()}-${server}`] = `webpack && grok publish ${server} --release`;
        }
        if (ts) Object.assign(_package.devDependencies, { 'ts-loader': 'latest', 'typescript': 'latest' });
        if (eslint) {
          Object.assign(_package.devDependencies, {
            'eslint': 'latest',
            'eslint-config-google': 'latest',
          }, ts ? {
            '@typescript-eslint/eslint-plugin': 'latest',
            '@typescript-eslint/parser': 'latest',
          } : {});
          Object.assign(_package.scripts, {
            'lint': `eslint src${ts ? ' --ext .ts' : ''}`,
            'lint-fix': `eslint src${ts ? ' --ext .ts' : ''} --fix`,
          });
        }
        if (jest) {
          Object.assign(_package.dependencies, {
            '@datagrok-libraries/utils': 'latest',
          });
          Object.assign(_package.scripts, {
            'test': 'grok test',
          });
        }

        // Save module names for installation prompt
        for (let [module, tag] of Object.entries(Object.assign({}, _package.dependencies, _package.devDependencies)))
          dependencies.push(`${module}@${tag}`);
        contents = JSON.stringify(_package, null, '\t');
      }
      if (file === 'package.js' && ts) copyFilePath = path.join(packageDir, 'package.ts');
      if (file === 'package-test.js' && ts) return false;
      if (file === 'package-test.ts' && !ts) return false;
      if (file === 'tsconfig.json' && !ts) return false;
      if (file === 'ts.webpack.config.js') return false;
      if (file === '.eslintrc.json') {
        if (!eslint) return false;
        if (ts) {
          let eslintConf = JSON.parse(contents);
          eslintConf.parser = '@typescript-eslint/parser';
          eslintConf.plugins = ['@typescript-eslint'];
          contents = JSON.stringify(eslintConf, null, '\t');
        }
      }
      if (file === 'gitignore') {
        copyFilePath = path.join(packageDir, '.gitignore');
        if (ts) contents += '\n# Emitted *.js files\nsrc/**/*.js\n';
      }
      if (file === 'npmignore') {
        copyFilePath = path.join(packageDir, '.npmignore');
      }
      fs.writeFileSync(copyFilePath, contents, 'utf8');
    } else if (stats.isDirectory()) {
      if (file === '.vscode' && !(ide == 'vscode' && platform == 'win32')) return;
      fs.mkdirSync(copyFilePath);
      // recursive call
      createDirectoryContents(name, config, origFilePath, copyFilePath, ide, ts, eslint, jest);
    }
  })
}

export function create(args: CreateArgs) {
  const options = ['ide', 'js', 'ts', 'eslint', 'jest'];
  const nOptions = Object.keys(args).length - 1;
  const nArgs = args['_'].length;

  if (nArgs > 2 || nOptions > 4) return false;
  if (nOptions && !Object.keys(args).slice(1).every(op => options.includes(op))) return false;
  if (args.js && args.ts) return color.error('Incompatible options: --js and --ts');
  const ts = !args.js && args.ts !== false;

  // Create `config.yaml` if it doesn't exist yet
  if (!fs.existsSync(grokDir)) fs.mkdirSync(grokDir);
  if (!fs.existsSync(confPath)) fs.writeFileSync(confPath, yaml.dump(confTemplate));

  const config = yaml.load(fs.readFileSync(confPath, { encoding: 'utf-8' })) as utils.Config;
  const confTest = validateConf(config);
  if (!confTest.value) {
    color.error(confTest.message);
    return false;
  }

  const name = nArgs === 2 ? args['_'][1] : curFolder;
  const validName = /^([A-Za-z\-_\d])+$/.test(name);

  if (validName) {
    let packageDir = curDir;
    let repositoryInfo: RepositoryInfo | null = null;

    if (curFolder !== name) {
      packageDir = path.join(packageDir, name);
      if (!fs.existsSync(packageDir)) {
        fs.mkdirSync(packageDir);
      }
    }
    if (!utils.isEmpty(packageDir)) {
      console.log();
      color.error('The package directory should be empty');
      return false;
    }

    exec('git rev-parse --is-inside-work-tree', { cwd: packageDir }, (err) => {
      if (err) return;
      const repository: RepositoryInfo = { type: 'git' };

      exec('git config --get remote.origin.url', { cwd: packageDir }, (err, stdout) => {
        if (err) return;
        repository.url = stdout.trim();

        exec('git rev-parse --show-prefix', { cwd: packageDir }, (err, stdout) => {
          if (err) return;
          const prefix = stdout.trim();
          repository.directory = prefix.endsWith('/') ? prefix.slice(0, -1) : prefix;

          if (repository.type && repository.url && repository.directory)
            repositoryInfo = repository;
        });
      });
    });

    process.on('beforeExit', () => {
      if (repositoryInfo) {
        const packagePath = path.join(packageDir, 'package.json');
        const p = JSON.parse(fs.readFileSync(packagePath, 'utf-8'));
        p.repository = repositoryInfo;
        fs.writeFileSync(packagePath, JSON.stringify(p, null, '\t'), 'utf-8');
      }
      process.exit();
    });

    createDirectoryContents(name, config, templateDir, packageDir, args.ide, ts, !!args.eslint, !!args.jest);
    color.success('Successfully created package ' + name);
    console.log(help.package(ts));
    console.log(`\nThe package has the following dependencies:\n${dependencies.join(' ')}\n`);
    console.log('Running `npm install` to get the required dependencies...\n');
    exec('npm install', { cwd: packageDir }, (err, stdout, stderr) => {
      if (err) throw err;
      else console.log(stderr, stdout);
    });
  } else {
    color.error('Package name may only include letters, numbers, underscores, or hyphens');
  }
  return true;
}

interface CreateArgs {
  _: string[],
  ide?: string,
  js?: boolean,
  ts?: boolean,
  eslint?: boolean,
  jest?: boolean
}

interface RepositoryInfo {
  type?: string,
  url?: string,
  directory?: string,
}
