import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
import { startDatagrok } from 'datagrok-api/datagrok';
import { readdir, writeFile } from 'fs/promises';
import { join, basename } from 'path';
declare let grok: typeof _grok, DG: typeof _DG;
export let _package: _DG.Package;

const testsExclude = [
    'sticky_meta.js',
    'benchmarks.js',
    'functions-annotations.js',
    'vector-functions-and-scripts.js'
];

async function main(): Promise<void> {
    const [apiUrl, token] = process.argv.slice(2);
    await startDatagrok({apiUrl: 'http://localhost:8888/api', apiToken: 'Bearer eyJhbGciOiJSUzI1NiIsInR5cCI6IkpXVCIsImtpZCI6ImF1dGgifQ.eyJpc3MiOiJodHRwOi8vbG9jYWxob3N0Ojg4ODgiLCJpZCI6IjA1OTlmNjUwLWFjNzAtMTFmMC1iMTNmLWQxMTZiN2U4YThjMiIsImV4cCI6IjIwMjUtMTEtMTFUMTM6MjI6MTUuMzUzODQxWiIsInN1YiI6eyJpZCI6Ijg3OGM0MmIwLTlhNTAtMTFlNi1jNTM3LTZiZjhlOWFiMDJlZSIsImxvZ2luIjoiYWRtaW4iLCJwcm9qZWN0Ijp7ImlkIjoiODc4YzQyYjAtOWE1MC0xMWU2LWM1MzctNmJmOGU5YWIwMjk5IiwibmFtZSI6IkFkbWluIn0sImdyb3VwIjp7ImlkIjoiYTRiNDU4NDAtOWE1MC0xMWU2LWM1MzctNmJmOGU5YWIwMmVlIn19fQ.e0sYhhuUNQGFSmB_ycqj_84RmwrgTyVOVWXrgMXat-CGjH6UlTzbR5rmSqi7SCa6wRQZd2Nul8BCKVraONrDd6d1teDUdwnFTALHt0AG7iq75jVzXVIqkML56VgHNhIMronmcv4KokQF-SHUa9J-Tiwl4BR1VqMN6vm85NQ1Rx6pR_8iW1KkKByBZPYUwxASJIFeM9vqexmIK9-hjtrj-VNF9BjBM_fyP00wqNh2S-gSMqPXKp-3ATYiuPM7N-SkR_-3CQUAW1_7YmoEkVermGFXlSIpQSGKiGoRMbJrwZ2nSkw6Pp4U8bUJMc2UkcTPDhvPYab1IhyoGKzOhYPo0g'});
    _package = await grok.dapi.packages.filter('shortName = "ApiTests"').first();
    if (!_package)
        throw new Error('ApiTests package should be installed');
    const dapiDir = join(__dirname, 'dapi');
    const files = await readdir(dapiDir);
    for (const file of files) {
        const baseName = basename(file);
        if (testsExclude.includes(baseName)) {
            console.log(`Skipping tests in ${baseName}`);
            continue;
        }
        if (baseName.endsWith('.js'))
           require(join(dapiDir, file));
    }
    const { runTests, tests} = await import('@datagrok-libraries/utils/src/test');
    const data = [];
    for (const key of Object.keys(tests))
        data.push(...(await runTests({
            category: key,
            test: undefined,
            testContext: undefined,
            nodeOptions: {package: _package}
        })))
    await writeFile('result.csv', DG.DataFrame.fromObjects(data)?.toCsv() ?? '', 'utf-8')
    console.log('Finished running tests');
    process.exit(0);
}


main().catch((err) => {
    console.error("Unhandled error:", err);
    process.exit(1);
});