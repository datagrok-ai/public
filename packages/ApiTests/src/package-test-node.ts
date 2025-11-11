import { readdir } from 'fs/promises';
import { writeFileSync } from 'fs';
import { basename, join } from 'path';

import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
import { startDatagrok } from 'datagrok-api/datagrok';

import yargs from 'yargs';
import { hideBin } from 'yargs/helpers';

declare let grok: typeof _grok, DG: typeof _DG;
export let _package: _DG.Package;

const testsExclude = [
    'sticky_meta.js',
    'benchmarks.js',
    'functions-annotations.js',
    'vector-functions-and-scripts.js'
];

async function main(): Promise<void> {
    const { apiUrl, devKey, concurrentRuns, categories, loop } = parseArgs();
    console.log('Exchanging devKey for token...');
    const apiToken = await getToken(apiUrl, devKey);
    console.log('Received token.');
    await startDatagrok({apiUrl, apiToken});
    _package = await grok.dapi.packages.filter('shortName = "ApiTests"').first();
    if (!_package)
        throw new Error('ApiTests package should be installed.');

    await loadTestFiles();

    const results: _DG.DataFrame[] = [];

    process.on('SIGINT', async () => {
        console.log('⚠️  Caught SIGINT (Ctrl+C)');
        await saveResults(results);
        process.exit(0);
    });

    process.on('SIGTERM', async () => {
        console.log('⚠️  Caught SIGTERM');
        await saveResults(results);
        process.exit(0);
    });

    await runPool(results, concurrentRuns, new Set(categories), loop);
    await saveResults(results);
    console.log('Finished running tests');
    process.exit(0);
}

function parseArgs() {
    const argv = yargs(hideBin(process.argv))
        .option('apiUrl', {
            type: 'string',
            describe: 'API URL of Datagrok server',
            demandOption: true,
        })
        .option('devKey', {
            type: 'string',
            describe: 'Developer key for authentication',
            demandOption: true,
        })
        .option('concurrentRuns', {
            type: 'number',
            describe: 'Number of concurrent runs to execute',
            default: 1,
        })
        .option('categories', {
            alias: 'c',
            type: 'string',
            array: true,
            describe: 'List of test categories to run. By default run all categories in dapi folder.',
            default: ['all'],
        })
        .option('loop', {
            alias: 'l',
            type: 'boolean',
            default: false,
            describe: 'If true, continuously repeat tasks with the specified number of concurrent runs',
        })
        .help()
        .parseSync();

    return {
        apiUrl: argv.apiUrl,
        devKey: argv.devKey,
        concurrentRuns: argv.concurrentRuns,
        categories: argv.categories,
        loop: argv.loop
    };
}

async function getToken(url: string, key: string) {
    const response = await fetch(`${url}/users/login/dev/${key}`, {method: 'POST'});
    const json = await response.json();
    if (json.isSuccess == true)
        return json.token;
    else
        throw new Error('Unable to login to server. Check your dev key.');
}

async function loadTestFiles(): Promise<void> {
    const dapiDir = join(__dirname, 'dapi');
    const files = await readdir(dapiDir);
    for (const file of files) {
        const baseName = basename(file);
        if (testsExclude.includes(baseName)) {
            console.log(`Skipping tests in ${baseName} because test file marked as unsupported.`);
            continue;
        }
        if (baseName.endsWith('.js'))
            require(join(dapiDir, file));
    }
}

async function run(concurrentRun: number, categories: Set<string>): Promise<_DG.DataFrame | null> {
    const { runTests, tests} = await import('@datagrok-libraries/utils/src/test');
    const data = [];
    const runAll = categories.has('all');
    for (const key of Object.keys(tests))
        if (runAll || categories.has(key))
            data.push(...(await runTests({
                category: key,
                test: undefined,
                testContext: undefined,
                nodeOptions: {package: _package}
            })));
    if (data.length === 0)
        return null;
    const df = DG.DataFrame.fromObjects(data)!;
    await df.columns.addNewCalculated('concurrent_run', `${concurrentRun}`, DG.COLUMN_TYPE.INT, false, false);
    return df;
}
const running = new Set();
async function runPool(results: _DG.DataFrame[], concurrency: number = 1, categories: Set<string>, infinite: boolean = false): Promise<void> {
    let counter = 0;

    function startNew() {
        const id = ++counter;
        const promise = run(id, categories)
            .then((r) => { if (r != null) results.push(r); })
            .catch((err: any) => console.error(`❌ Task ${id} failed:`, err))
            .finally(() => {
                running.delete(promise);
                if (infinite)
                    startNew();
            });
        running.add(promise);
    }

    for (let i = 0; i < concurrency; i++)
        startNew();

    if (infinite)
        await new Promise(() => {});
    else
        await Promise.allSettled(running);
}

async function saveResults(results: _DG.DataFrame[]): Promise<void> {
    console.log('Creating test report...')
    if (results.length === 0) {
        console.log('Result are empty');
        return;
    }
    const res: _DG.DataFrame = results[0];
    for (let i = 1; i < results.length; i++)
        res.append(results[i], true);
    await res.columns.addNewCalculated('stress_test', 'true', DG.COLUMN_TYPE.BOOL, false, false);
    await res.columns.addNewCalculated('benchmark', 'false', DG.COLUMN_TYPE.BOOL, false, false);
    console.log('Writing results in test-report.csv.')
    writeFileSync('test-report.csv', res.toCsv(), 'utf-8')
    console.log('Saved results in test-report.csv.')
}

main().catch((err) => {
    console.error("Unhandled error:", err);
    process.exit(1);
});
