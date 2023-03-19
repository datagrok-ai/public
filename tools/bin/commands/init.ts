import fs from 'fs';
import os from 'os';
import path from 'path';
import yaml from 'js-yaml';
import * as utils from '../utils/utils';
import * as color from '../utils/color-utils';
import { validateConf } from '../validators/config-validator';


const curDir = process.cwd();
const platform = os.platform();
const grokDir = path.join(os.homedir(), '.grok');
const confPath = path.join(grokDir, 'config.yaml');
const packageJsonPath = path.join(curDir, 'package.json');
const tsConfigPath = path.join(curDir, 'tsconfig.json');
const webpackConfigPath = path.join(curDir, 'webpack.config.js');
const templateDir = path.join(path.dirname(path.dirname(__dirname)), 'package-template');


export function init(args: InitArgs) {
  const options = ['ide', 'eslint', 'test', 'ts'];
  const nOptions = Object.keys(args).length - 1;
  const nArgs = args['_'].length;

  if (nArgs > 1 || nOptions > options.length) {
    color.error('Incorrect number of arguments and options.');
    return false;
  }

  const passedOptions = Object.keys(args).slice(1);
  if (nOptions) {
    let hasUnknownOpt = false;
    for (const op of passedOptions) {
      if (!options.includes(op)) {
        if (op !== 'h' && op !== 'help')
          color.error(`Unknown option: ${op}`);
        hasUnknownOpt = true;
      }
    }
    if (hasUnknownOpt)
      return false;
  }

  if (!utils.isPackageDir(curDir)) {
    color.error('File `package.json` not found. Run the command from the package directory');
    return false;
  }

  if (args.ts) {
    fs.writeFileSync(webpackConfigPath, fs.readFileSync(path.join(templateDir, 'ts.webpack.config.js')));
    fs.writeFileSync(tsConfigPath, fs.readFileSync(path.join(templateDir, 'tsconfig.json')));
    const packageJsPath = path.join(curDir, 'src', 'package.js');
    fs.writeFileSync(path.join(curDir, 'src', 'package.ts'), fs.existsSync(packageJsPath) ?
      fs.readFileSync(packageJsPath) : fs.readFileSync(path.join(templateDir, 'src', 'package.ts')));
    const packageData = JSON.parse(fs.readFileSync(packageJsonPath, 'utf-8'));
    Object.assign(packageData.devDependencies, { 'ts-loader': 'latest', 'typescript': 'latest' });
    if ('eslint' in packageData.devDependencies)
      Object.assign(packageData.devDependencies, {
        '@typescript-eslint/eslint-plugin': 'latest',
        '@typescript-eslint/parser': 'latest',
      });
    fs.writeFileSync(packageJsonPath, JSON.stringify(packageData, null, 2), 'utf-8');
    color.success('TypeScript support has been added.');
  }

  if (args.ide) {
    const config = yaml.load(fs.readFileSync(confPath, { encoding: 'utf-8' })) as utils.Config;
    const confTest = validateConf(config);
    if (!confTest.value) {
      color.warn('Invalid configuration. Skipping `ide`...');
      color.error(confTest.message);
    }

    if (args.ide === 'vscode') {
      const ideConfPath = path.join(curDir, '.vscode');
      const templateConfPath = path.join(templateDir, '.vscode');
      if (!fs.existsSync(ideConfPath))
        fs.mkdirSync(ideConfPath);
      const files = fs.readdirSync(templateConfPath);
      for (const file of files) {
        let contents = fs.readFileSync(path.join(templateConfPath, file), 'utf-8');
        if (file === 'tasks.json' && platform !== 'win32')
          contents = contents.replace(/(?<="command": ").*(?=")/, 'webpack && grok publish #{GROK_HOST_ALIAS}');
        contents = contents.replace(/#{GROK_HOST_ALIAS}/g, config.default);
        contents = contents.replace(/#{GROK_HOST}/g, /localhost|127\.0\.0\.1/.test(
          config['servers'][config.default]['url']) ? 'http://localhost:63343/login.html'
          : (new URL(config['servers'][config.default]['url'])).origin);
        fs.writeFileSync(path.join(ideConfPath, file), contents, 'utf-8');
        color.success('IDE configuration has been added.');
      }
    }
  }

  if (args.eslint) {
    const packageData = JSON.parse(fs.readFileSync(packageJsonPath, 'utf-8'));
    Object.assign(packageData.devDependencies, {
      'eslint': 'latest',
      'eslint-config-google': 'latest',
    });
    const ts = fs.existsSync(tsConfigPath);
    if (ts) {
      Object.assign(packageData.devDependencies, {
        '@typescript-eslint/eslint-plugin': 'latest',
        '@typescript-eslint/parser': 'latest',
      });
    }
    Object.assign(packageData.scripts, {
      'lint': `eslint src${ts ? ' --ext .ts' : ''}`,
      'lint-fix': `eslint src${ts ? ' --ext .ts' : ''} --fix`,
    });
    fs.writeFileSync(packageJsonPath, JSON.stringify(packageData, null, 2), 'utf-8');
    const eslintConfigPath = path.join(curDir, '.eslintrc.json');
    const eslintTemplatePath = path.join(templateDir, '.eslintrc.json');
    if (ts) {
      let contents = fs.readFileSync(eslintTemplatePath, 'utf-8');
      const eslintConf = JSON.parse(contents);
      eslintConf.parser = '@typescript-eslint/parser';
      eslintConf.plugins = ['@typescript-eslint'];
      contents = JSON.stringify(eslintConf, null, 2);
      fs.writeFileSync(eslintConfigPath, contents, 'utf-8');
    } else
      fs.writeFileSync(eslintConfigPath, eslintTemplatePath);
    color.success('Linter configuration has been added.');
  }

  if (args.test) {
    const tsPath = path.join(curDir, 'src', 'package.ts');

    if (!fs.existsSync(tsPath)) {
      color.warn('Tests can only be added to TypeScript packages. Skipping `test`...');
    } else if (!fs.existsSync(webpackConfigPath)) {
      color.warn('Webpack configuration not found. Skipping `test`...');
    } else {
      const config = fs.readFileSync(webpackConfigPath, 'utf-8');
      if (!/(?<=entry:\s*{\s*(\r\n|\r|\n))[^}]*test:/.test(config)) {
        const entryIdx = config.search(/(?<=entry:\s*{\s*(\r\n|\r|\n)).*/);
        if (entryIdx === -1)
          color.error('Entry point not found during webpack config parsing');
        else {
          const testEntry = "    test: {filename: 'package-test.js', library: " +
            "{type: 'var', name:`${packageName}_test`}, import: './src/package-test.ts'},";
          fs.writeFileSync(webpackConfigPath, config.slice(0, entryIdx) + testEntry +
            config.slice(entryIdx), 'utf-8');
        }
      }

      const packageData = JSON.parse(fs.readFileSync(packageJsonPath, 'utf-8'));
      Object.assign(packageData.dependencies, {
        '@datagrok-libraries/utils': 'latest',
      });
      Object.assign(packageData.scripts, {
        'test': 'grok test',
      });
      fs.writeFileSync(packageJsonPath, JSON.stringify(packageData, null, 2), 'utf-8');

      const packageTestPath = path.join(curDir, 'src', 'package-test.ts');
      if (!fs.existsSync(packageTestPath))
        fs.writeFileSync(packageTestPath, fs.readFileSync(path.join(templateDir, 'src', 'package-test.ts')));

      const testsDir = path.join(curDir, 'src', 'tests');
      if (!fs.existsSync(testsDir))
        fs.mkdirSync(testsDir);

      fs.writeFileSync(path.join(testsDir, 'test-examples.ts'),
        fs.readFileSync(path.join(path.dirname(templateDir), 'entity-template', 'test.ts')));

      fs.writeFileSync(packageTestPath, `import './tests/test-examples';\n` +
        fs.readFileSync(packageTestPath, 'utf-8'), 'utf-8');
      color.success('Tests support has been added.');
    }
  }

  return true;
}

interface InitArgs {
  _: string[],
  ide?: string,
  eslint?: boolean,
  test?: boolean,
  ts?: boolean,
}
