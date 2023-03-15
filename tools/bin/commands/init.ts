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
      }
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
