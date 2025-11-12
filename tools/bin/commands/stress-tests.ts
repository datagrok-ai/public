import fs from 'fs';
import path from 'path';
import * as utils from '../utils/utils';
import * as testUtils from '../utils/test-utils';
import {getDevKey} from "../utils/test-utils";

const cwd = process.cwd();

export async function stressTests(args: StressTestArgs): Promise<boolean> {
    const pkgPath = path.join(cwd, 'package.json');
    if (!fs.existsSync(pkgPath)) {
        console.error('❌ Error: This command must be executed from the ApiTests package folder (package.json not found).');
        process.exit(1);
    }
    let pkg;
    try {
        pkg = JSON.parse(fs.readFileSync(pkgPath, 'utf8'));
    } catch (e) {
        console.error('❌ Error: Failed to read package.json. Make sure it is valid JSON.');
        process.exit(1);
    }
    if (pkg.name !== '@datagrok/api-tests') {
        console.error(`❌ Error: This command must be executed from the ApiTests package folder. Found package name: '${pkg.name ?? 'undefined'}'`);
        process.exit(1);
    }
    const config = getDevKey(args.host);
    await testUtils.loadPackage('', 'ApiTests', args.host, args['skip-publish'], args['skip-build']);
    await utils.runScript(`npm run build-node`, '');
    try {
        await utils.runScript(`node -r ./tsconfig-paths-bootstrap.js dist-node/package-test-node.js --apiUrl=${config.url} --devKey=${config.key} --concurrentRuns=${args["concurrent-runs"] ?? 1}${args.loop ? ' --loop' : ''}`, '');
        return true;
    } catch (e) {
        console.error(`❌ Error: Something went wrong: ${e}`);
        return false;
    }
}

interface StressTestArgs {
    _: string[];
    host: string;
    'skip-build'?: boolean;
    'skip-publish'?: boolean;
    'concurrent-runs'?: number;
    loop?: number;
}