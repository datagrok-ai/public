import { readdir } from 'fs/promises';
import { writeFileSync } from 'fs';
import { fork, ChildProcess } from 'child_process';
import { availableParallelism } from 'os';
import { basename, dirname, join } from 'path';
import { fileURLToPath, pathToFileURL } from 'url';

import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';

import yargs from 'yargs';
import { hideBin } from 'yargs/helpers';

import { setTestPackage } from './test-package';

declare let grok: typeof _grok, DG: typeof _DG;
export let _package: _DG.Package;

// Entries are basenames or dir-qualified paths (e.g. 'shell/events.ts')
const testsExclude = [
    'sticky_meta.ts',
    'benchmarks.ts',
    'functions-annotations.ts',
    'vector-functions-and-scripts.ts',
    'cache.ts',          // client function-result cache: IndexedDB + grok.functions.register (browser-only)
    'shell/events.ts',   // shell workspace events (table/view add) — browser-only
    'shell/settings.ts', // client shell settings object — browser-only
    'shell/ml.ts',       // registers client JS funcs at import (grok.functions.register)
    'dataframe/add-new-column.ts', // dialog-driven (PowerPack addNewColumnDialog)
    'dataframe/calculated-columns-bench.ts', // table views and layouts
    'dapi/domain-frame-editor.ts', // DG.DomainFrameEditor + DG.Grid.attachEditor (browser-only)
];

// UI-independent test folders that run under Node; grid/, widgets/, ai/, packages/
// and utils/ stay browser-only (Dart client / DOM dependencies).
const nodeTestDirs = ['dapi', 'dataframe', 'functions', 'bitset', 'valuematcher', 'property', 'stats', 'shell'];

// Sibling packages whose tests join this suite, so the sweep covers more than
// the platform API — DBTests contributes read-only Postgres queries through grok_connect.
// Each package resolves its own copy of the test library, so their categories land in a
// different registry object and have to be merged into the one this runner filters.
const extraTestPackages = ['DBTests'];

// 'functional': run every loaded test (UI-dependent ones self-skip via skipReason);
// 'stress': only stressTest-marked tests — the concurrency-sweep baseline.
// Defaults to 'stress': `grok stresstest` (the nightly Stress-Tests baseline) invokes
// this runner without a mode flag, and the curated stress set must not silently widen.
// Functional runs opt in explicitly (npm run start-node passes --mode=functional).
let mode: 'functional' | 'stress' = 'stress';

async function main(): Promise<void> {
    const { apiUrl, devKey, token, concurrentRuns, categories, loop, concurrencyRange, workers, worker } = parseArgs();
    const apiToken = token ?? await getToken(apiUrl, devKey);
    console.log('Received token.');
    // Loaded dynamically (not a static import) so Node/tsx resolve datagrok-api's
    // named exports through cjs-module-lexer instead of failing at link time.
    const { startDatagrok } = await import('datagrok-api/datagrok');
    await startDatagrok({apiUrl, apiToken});
    // Bind datagrok-api/{dg,grok,ui} to the runtime objects startDatagrok() just produced,
    // so every test file shares one copy of the js-api classes. Must precede any test-file
    // import. Same module instance the --import already loaded.
    const loader = await import(pathToFileURL(
        join(dirname(fileURLToPath(import.meta.url)), '..', 'node-test-loader', 'register.mjs')).href);
    loader.bindRuntimeGlobals();
    _package = await grok.dapi.packages.filter('shortName = "ApiTests"').first();
    if (!_package)
        throw new Error('ApiTests package should be installed.');
    setTestPackage(_package);

    await loadTestFiles();

    if (worker)
        return await workerMain(new Set(categories));

    const results: _DG.DataFrame[] = [];
    const pool = workers > 1 ? await WorkerPool.start(workers, results) : null;

    process.on('SIGINT', async () => {
        console.log('⚠️  Caught SIGINT (Ctrl+C)');
        await pool?.stop();
        await saveResults(results);
        process.exit(0);
    });

    process.on('SIGTERM', async () => {
        console.log('⚠️  Caught SIGTERM');
        await pool?.stop();
        await saveResults(results);
        process.exit(0);
    });

    process.on("uncaughtException", (err) => {
        console.error("Uncaught exception:", err);
    });

    process.on("unhandledRejection", (reason) => {
        console.error("Unhandled rejection:", reason);
    });

    if (concurrencyRange) {
        for (const runs of concurrencyRange)
            await (pool ? pool.runLevel(runs) : runPool(results, runs, new Set(categories)));
    }
    else
        await (pool ? pool.runLevel(concurrentRuns, loop) : runPool(results, concurrentRuns, new Set(categories), loop));

    await pool?.stop();
    await saveResults(results);
    console.log('Finished running tests');
    if (loadFailures.length > 0) {
        console.error(`\n❌ ${loadFailures.length} test file(s) failed to load: ${loadFailures.join(', ')}`);
        process.exit(1);
    }
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
            describe: 'Developer key for authentication (deprecated - prefer --token)',
        })
        .option('token', {
            type: 'string',
            describe: 'Session token, e.g. from `grok s token`. Preferred over --devKey',
        })
        .option('concurrentRuns', {
            type: 'number',
            describe: 'Number of concurrent runs to execute'
        })
        .coerce('concurrencyRange', (value) => {
            if (!value) return undefined;

            const match = value.trim().match(/^(\d+)\s*-\s*(\d+)$/);
            if (!match)
                throw new Error(`Invalid format "${value}". Expected format: "start-end", e.g. "1-10".`);

            const start = Number(match[1]);
            const end = Number(match[2]);

            if (start <= 0 || end <= 0)
                throw new Error(`Range must contain positive integers.`);
            if (start >= end)
                throw new Error(`Invalid range "${value}". Start must be < end.`);

            return { start, end };
        })
        .option('categories', {
            alias: 'c',
            type: 'string',
            array: true,
            describe: 'List of test categories to run. By default run all categories in dapi folder.',
            default: ['all'],
        })
        .option('workers', {
            type: 'number',
            describe: 'Processes the concurrent runs are spread over. Default: one per 16 runs of the widest ' +
                'level, at most a quarter of the cores; 1 keeps every run in this process'
        })
        .option('worker', {
            type: 'boolean',
            hidden: true,
            describe: 'Internal: this process runs the share of the runs its parent sends it'
        })
        .option('loop', {
            alias: 'l',
            type: 'boolean',
            default: false,
            describe: 'If true, continuously repeat tasks with the specified number of concurrent runs',
        })
        .option('mode', {
            type: 'string',
            choices: ['functional', 'stress'] as const,
            default: 'stress',
            describe: 'functional: run all loaded tests once; stress: only stressTest-marked tests (concurrency baseline)',
        })
        .option('step', {
            type: 'number',
            describe: 'Step to use for concurrencyRange'
        })
        .conflicts('concurrencyRange', ['concurrentRuns'])
        .conflicts('concurrentRuns', 'concurrencyRange')
        .check((argv) => {
            if (argv.loop === true && argv.concurrencyRange !== undefined)
                throw new Error('You cannot use --loop together with --concurrencyRange.');
            if (argv.step && argv.concurrencyRange === undefined)
                throw new Error('The --step parameter can only be used with --concurrencyRange.');

            if (argv.step && argv.step <= 0)
                throw new Error('The --step parameter must be a positive integer.');

            return true;
        })
        .help()
        .parseSync();

    mode = argv.mode as typeof mode;
    const res: any = {
        apiUrl: argv.apiUrl,
        devKey: argv.devKey,
        token: argv.token,
        concurrentRuns: argv.concurrentRuns,
        categories: argv.categories,
        loop: argv.loop
    };
    if (argv.concurrencyRange) {
        const { start, end } = argv.concurrencyRange;
        const step = argv.step ?? 1;

        const values = [];
        for (let i = start; i <= end; i += step)
            values.push(i);

        res.concurrencyRange = values;
    }

    if (!res.concurrentRuns && !res.concurrencyRange)
        res.concurrentRuns = 1;
    // One process runs every test of every run on one thread: past a few dozen runs the tests wait for
    // that thread, not for the server, and the load test measures its own runner.
    const widest = res.concurrencyRange ? Math.max(...res.concurrencyRange) : res.concurrentRuns;
    res.worker = argv.worker === true;
    res.workers = Math.max(1, Math.floor(argv.workers ??
        Math.min(Math.ceil(widest / 16), Math.floor(availableParallelism() / 4))));

    return res;
}

async function getToken(url: string, key: string) {
    if (!key)
        throw new Error('No credentials: pass --token (from `grok s token`) or --devKey');
    console.log('Exchanging devKey for token...');
    // Raw fetch is intentional here: this runs before startDatagrok(), so the grok
    // client (and grok.dapi) isn't initialized yet — there's no dapi layer to use.
    const response = await fetch(`${url}/users/login/dev`, {method: 'POST', headers: {'Authorization': `Dev ${key}`}});
    const json = await response.json();
    if (json.isSuccess == true)
        return json.token;
    else
        throw new Error('Unable to login to server. Check your dev key.');
}

const loadFailures: string[] = [];

async function loadTestFiles(): Promise<void> {
    const srcDir = dirname(fileURLToPath(import.meta.url));
    for (const dir of nodeTestDirs) {
        const testDir = join(srcDir, dir);
        const files = await readdir(testDir);
        for (const file of files) {
            const baseName = basename(file);
            if (testsExclude.includes(baseName) || testsExclude.includes(`${dir}/${baseName}`)) {
                console.log(`Skipping tests in ${dir}/${baseName} because test file marked as unsupported.`);
                continue;
            }
            if (baseName.endsWith('.ts')) {
                try {
                    await import(pathToFileURL(join(testDir, file)).href);
                } catch (e: any) {
                    loadFailures.push(`${dir}/${baseName}`);
                    console.error(`❌ Failed to load ${dir}/${baseName} (its tests will not run): ${e?.message ?? e}`);
                }
            }
        }
    }
    const {tests} = await import('@datagrok-libraries/test/src/test');
    for (const pkg of extraTestPackages) {
        try {
            const entry = pathToFileURL(join(srcDir, '..', '..', pkg, 'src', 'package-test.ts')).href;
            const loaded: any = await import(entry);
            const added = Object.keys(loaded.tests ?? {});
            Object.assign(tests, loaded.tests);
            console.log(`Loaded ${added.length} categories from ${pkg}: ${added.join(', ')}`);
        } catch (e: any) {
            loadFailures.push(pkg);
            console.error(`❌ Failed to load ${pkg} (its tests will not run): ${e?.message ?? e}`);
        }
    }
}

async function run(concurrentRun: number, categories: Set<string>, concurrency: number): Promise<_DG.DataFrame | null> {
    const data = await runRows(categories);
    return data == null ? null : await rowsToFrame(data, concurrentRun, concurrency);
}

/** The result rows of one run of the suite; null when no test ran. */
async function runRows(categories: Set<string>): Promise<any[] | null> {
    const { runTests, tests} = await import('@datagrok-libraries/test/src/test');
    const data = [];
    const runAll = categories.has('all');
    // The stress suite is the stressTest-marked tests. Filter the registry directly
    // (not all installed @datagrok-libraries/test versions honor the stressTest run
    // option), so only stress-marked tests run and the baseline can be 100%.
    // Functional mode runs everything that was loaded.
    if (mode === 'stress') {
        for (const key of Object.keys(tests)) {
            const cat: any = (tests as any)[key];
            if (cat?.tests)
                cat.tests = cat.tests.filter((t: any) => t.options?.stressTest);
        }
    }
    // Wall-clock window for this task, so stress runs can be correlated with the
    // container CPU/memory time series (UTC, matching the docker_stats collector).
    const taskStart = new Date().toISOString();
    for (const key of Object.keys(tests))
        if ((runAll || categories.has(key)) && (tests as any)[key]?.tests?.length)
            data.push(...(await runTests({
                category: key,
                test: undefined,
                testContext: undefined,
                stressTest: mode === 'stress',
                nodeOptions: {package: _package}
            })));
    if (data.length === 0)
        return null;
    const taskFinish = new Date().toISOString();
    for (const row of data as any[]) {
        row.task_start = taskStart;
        row.task_finish = taskFinish;
    }
    return data;
}

async function rowsToFrame(data: any[], concurrentRun: number, concurrency: number): Promise<_DG.DataFrame> {
    const df = DG.DataFrame.fromObjects(data)!;
    await df.columns.addNewCalculated('concurrent_run', `${concurrentRun}`, DG.COLUMN_TYPE.INT, false, false);
    await df.columns.addNewCalculated('total_concurrent_runs', `${concurrency}`, DG.COLUMN_TYPE.INT, false, false);
    return df;
}
const running = new Set();
async function runPool(results: _DG.DataFrame[], concurrency: number = 1, categories: Set<string>, infinite: boolean = false): Promise<void> {
    let counter = 0;

    function startNew() {
        const id = ++counter;
        const promise = run(id, categories, concurrency)
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

/** In a worker process: runs the share of each level its parent sends and hands the rows of every run back. */
async function workerMain(categories: Set<string>): Promise<void> {
    process.on('SIGINT', () => {});
    process.on('SIGTERM', () => {});
    process.on('disconnect', () => process.exit(0));
    process.on('message', (m: any) => {
        if (m.type === 'stop')
            process.exit(0);
        if (m.type !== 'level')
            return;
        let left = m.ids.length;
        let restarts = 0;
        const start = (id: number) => {
            runRows(categories)
                .then((rows) => { if (rows != null) process.send!({type: 'rows', run: id, total: m.total, rows}); })
                .catch((err: any) => console.error(`❌ Task ${id} failed:`, err))
                .finally(() => {
                    if (m.infinite)
                        start(m.total + 1 + m.worker + m.workers * restarts++);
                    else if (--left === 0)
                        process.send!({type: 'done'});
                });
        };
        for (const id of m.ids)
            start(id);
    });
    process.send!({type: 'ready'});
    await new Promise(() => {});
}

/** The processes the runs of a level are spread over, started once and used for every level. */
class WorkerPool {
    private frames: Promise<void> = Promise.resolve();
    private levelDone: (() => void)[] = [];

    private constructor(private children: ChildProcess[], private results: _DG.DataFrame[]) {}

    static async start(count: number, results: _DG.DataFrame[]): Promise<WorkerPool> {
        console.log(`Starting ${count} worker processes...`);
        const children = Array.from({length: count}, () => fork(fileURLToPath(import.meta.url),
            [...process.argv.slice(2), '--worker'], {execArgv: process.execArgv, serialization: 'advanced'}));
        const pool = new WorkerPool(children, results);
        await Promise.all(children.map((child, i) => new Promise<void>((resolve, reject) => {
            child.on('message', (m: any) => {
                if (m.type === 'ready')
                    resolve();
                else if (m.type === 'rows')
                    pool.frames = pool.frames.then(async () => { results.push(await rowsToFrame(m.rows, m.run, m.total)); })
                        .catch((e: any) => console.error(`❌ The results of run ${m.run} were lost:`, e));
                else if (m.type === 'done')
                    pool.levelDone[i]?.();
            });
            child.on('exit', (code) => {
                if (code !== 0 && code !== null)
                    console.error(`❌ Worker ${i} exited with code ${code}`);
                pool.levelDone[i]?.();
                reject(new Error(`Worker ${i} exited before it was ready`));
            });
        })));
        console.log('Workers are ready.');
        return pool;
    }

    /** Runs [runs] concurrent runs, as evenly over the workers as they divide; an infinite level never ends. */
    async runLevel(runs: number, infinite: boolean = false): Promise<void> {
        const count = this.children.length;
        const waits: Promise<void>[] = [];
        let next = 1;
        for (let i = 0; i < count; i++) {
            const share = Math.floor(runs / count) + (i < runs % count ? 1 : 0);
            if (share === 0)
                continue;
            const ids = Array.from({length: share}, () => next++);
            waits.push(new Promise<void>((resolve) => this.levelDone[i] = resolve));
            this.children[i].send({type: 'level', ids, total: runs, infinite, worker: i, workers: count});
        }
        await Promise.all(waits);
        await this.frames;
    }

    /** Stops the workers; the rows they already sent are in the results when this returns. */
    async stop(): Promise<void> {
        for (const child of this.children)
            if (child.connected)
                child.send({type: 'stop'});
        await this.frames;
    }
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
    await res.columns.addNewCalculated('stress_test', `${mode === 'stress'}`, DG.COLUMN_TYPE.BOOL, false, false);
    await res.columns.addNewCalculated('benchmark', 'false', DG.COLUMN_TYPE.BOOL, false, false);
    console.log('Writing results in test-report.csv.')
    writeFileSync('test-report.csv', res.toCsv(), 'utf-8')
    console.log('Saved results in test-report.csv.')
}

main().catch((err) => {
    console.error("Unhandled error:", err);
    process.exit(1);
});
