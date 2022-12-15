import { exec } from 'child_process';
import fs from 'fs';
import os from 'os';
import path from 'path';
import yaml from 'js-yaml';
import * as utils from '../utils/utils';
import * as color from '../utils/color-utils';


export function test(args: TestArgs): boolean {
  const options = Object.keys(args).slice(1);
  const nArgs = args['_'].length;
  const curDir = process.cwd();
  const grokDir = path.join(os.homedir(), '.grok');
  const confPath = path.join(grokDir, 'config.yaml');

  if (nArgs > 1 || options.length > 1 || (options.length === 1 && options[0] !== 'host'))
    return false;

  if (!utils.isPackageDir(curDir)) {
    color.error('File `package.json` not found. Run the command from the package directory');
    return false;
  }

  if (!fs.existsSync(confPath)) {
    color.error(`File \`${confPath}\` not found. Run \`grok config\` to set up the config file`);
    return false;
  }

  const config = yaml.load(fs.readFileSync(confPath, { encoding: 'utf-8' })) as utils.Config;

  if (args.host) {
    if (args.host in config.servers) {
      process.env.HOST = args.host;
      console.log('Environment variable `HOST` is set to', args.host);
    } else {
      color.error(`Unknown server alias. Please add it to ${confPath}`);
      return false;
    }
  } else if (config.default) {
    process.env.HOST = config.default;
    console.log('Environment variable `HOST` is set to', config.default);
  }

  const packageData = JSON.parse(fs.readFileSync(path.join(curDir, 'package.json'), { encoding: 'utf-8' }));
  const fullName = packageData.friendlyName || packageData.fullName;
  if (fullName) {
    process.env.TARGET_PACKAGE = fullName;
    console.log('Environment variable `TARGET_PACKAGE` is set to', process.env.TARGET_PACKAGE);
  } else {
    color.error('Invalid full package name. Set `friendlyName` or `fullName` field in `package.json`');
    return false;
  }

  color.info(`Publishing package "${process.env.TARGET_PACKAGE}" to ${process.env.HOST}...`);
  exec(`grok publish ${process.platform === 'win32' ? '%HOST%' : '${HOST}'}`, (err, stdout, stderr) => {
    if (err)
      throw err;
    else {
      console.log(stdout);
      color.error(stderr);
    }

    color.info('Starting tests...');
    exec('npm run test', {cwd: path.dirname(path.dirname(__dirname))}, (err, stdout, stderr) => {
      if (err) throw err;
      else {
        console.log(stdout);
        console.log(stderr);
      }
    });
  });

  return true;
}

interface TestArgs {
  _: string[],
  host?: string,
}
