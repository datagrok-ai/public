import { runWorker } from "../utils/test-utils";
import * as utils from '../utils/utils';
import fs from 'fs';
import path from 'path';
import os from 'os';
import yaml from 'js-yaml';

const grokDir = path.join(os.homedir(), '.grok');
const confPath = path.join(grokDir, 'config.yaml');

export  async function loadNPM(args: LoadNPMArgs) {
    const config = yaml.load(fs.readFileSync(confPath, { encoding: 'utf-8' })) as utils.Config;
    utils.setHost(args.host, config);
    await runWorker([{
        package: 'Loading Npm Packages',
        params: {
            test: '',
            category: '',
            options: {
                catchUnhandled: undefined,
                report: undefined
            }
        }
    }], {}, 0);
    process.exit(0);
}

interface LoadNPMArgs {
    host?: string;
}