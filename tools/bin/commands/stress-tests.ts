import fs from 'fs';
import { spawn } from 'child_process';
import path from 'path';
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
    await testUtils.loadPackage('', 'ApiTests', args.host, args['skip-publish'], args['skip-build'], false, true);
    try {
        await run(config, args);
        return true;
    } catch (e) {
        console.error(`❌ Error: Something went wrong: ${e}`);
        return false;
    }
}

async function run(config: { url: string, key: string }, args: StressTestArgs): Promise<number> {
    const processArgs: string[] = [];
    // Run the TypeScript runner directly via tsx (no build step); the asset/dayjs shims
    // in node-test-loader let the browser-oriented sources load under Node.
    processArgs.push('--import');
    processArgs.push('tsx');
    processArgs.push('--import');
    processArgs.push('./node-test-loader/register.mjs');
    processArgs.push('src/package-test-node.ts');
    processArgs.push(`--apiUrl=${config.url}`);
    processArgs.push(`--devKey=${config.key}`);
    if (args["concurrent-runs"])
        processArgs.push(`--concurrentRuns=${args["concurrent-runs"]}`);
    if (args['loop'])
        processArgs.push(`--loop`);
    if (args['concurrency-range'])
        processArgs.push(`--concurrencyRange=${args["concurrency-range"]}`);
    if (args['step'])
        processArgs.push(`--step=${args["step"]}`);

    const child = spawn(
        'node',
        processArgs,
        {
            cwd: '',
            stdio: ['inherit', 'inherit', 'inherit'],
        }
    );

    return new Promise((resolve, reject) => {
        child.on('exit', (code) => {
            if (code === 0) resolve(code);
            else reject(new Error(`Stress tests exited with code ${code}`));
        });
    });
}

interface StressTestArgs {
    _: string[];
    host: string;
    'skip-build'?: boolean;
    'skip-publish'?: boolean;
    'concurrent-runs'?: number;
    loop?: number;
    'concurrency-range'?: string;
    step?: number;
}