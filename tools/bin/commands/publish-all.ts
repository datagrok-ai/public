import { loadPackages } from "../utils/test-utils";
import * as utils from '../utils/utils';
import yaml from 'js-yaml';
import fs from 'fs';
import os from 'os';
import path from 'path';

const { exec } = require('child_process');

const curDir = process.cwd();
const grokDir = path.join(os.homedir(), '.grok');
const confPath = path.join(grokDir, 'config.yaml');


export async function publishAll(args: PublishAllArgs) {
    const config = yaml.load(fs.readFileSync(confPath, { encoding: 'utf-8' })) as utils.Config;

    utils.setHost(args.host, config);;
    await loadPackages(curDir, 'all', args.host, false, args["skip-build"], args["link-package"], args.release);
}

interface PublishAllArgs {
  _: string[],
  host?: string,
  release?: boolean, 
  'link-package'?: boolean;
  'skip-build'?: boolean;
}
