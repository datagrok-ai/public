import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
import {startDatagrok} from 'datagrok-api/datagrok';
import fs from "fs";
declare let grok: typeof _grok, DG: typeof _DG;

export let _package: _DG.Package;

const testsExclude = [
    'benchmarks.js'
];

async function main(): Promise<void> {
    const [apiUrl, token] = process.argv.slice(2);
    await startDatagrok({apiUrl: 'http://localhost:8888/api', apiToken: 'Bearer eyJhbGciOiJSUzI1NiIsInR5cCI6IkpXVCIsImtpZCI6ImF1dGgifQ.eyJpc3MiOiJodHRwOi8vbG9jYWxob3N0Ojg4ODgiLCJpZCI6ImUyZDJhMGYwLWFjNzAtMTFmMC1iMjZkLTMzYjQyYTVhY2VkZSIsImV4cCI6IjIwMjUtMTEtMTBUMTA6NDc6NDQuNzE2ODE0WiIsInN1YiI6eyJpZCI6Ijg3OGM0MmIwLTlhNTAtMTFlNi1jNTM3LTZiZjhlOWFiMDJlZSIsImxvZ2luIjoiYWRtaW4iLCJwcm9qZWN0Ijp7ImlkIjoiODc4YzQyYjAtOWE1MC0xMWU2LWM1MzctNmJmOGU5YWIwMjk5IiwibmFtZSI6IkFkbWluIn0sImdyb3VwIjp7ImlkIjoiYTRiNDU4NDAtOWE1MC0xMWU2LWM1MzctNmJmOGU5YWIwMmVlIn19fQ.fulNJbDWTTTzAITXtVwOlEfrOpZATRa5Ih5yHB2j-U0CNe8QNAFDlX3s95W-0jX2PM9PE8sXnQQZrbcUufjSj31J6Fcj-Mmrtn_eHqIX6ShTMl79LIR1hmBobgSpfmAMx2vXX-YjeSLqyyZ8KB4184nJKxW9w9TGF660mG8U_IxqHHQqbACqc2WjIYD0kJGV_zGYZkE92H4HXviChSNGxPn4W7LaCBFI958egd5qjMuvTtC2IEg0lZ1Dnw5UROcZJQ0rvh8jZeAx8AWQ0scKYtfPwjnKAyuZiNNblYFfrvyUfoP0cF67D8ch-p3mhbdFCkUyYk2HpeP3CzBB-XtO2Q'});
    _package = await grok.dapi.packages.filter('shortName = "ApiTests"').first();
    if (!_package)
        throw new Error('ApiTests package should be installed');
    const { readdir, writeFile } = await import('fs/promises')
    const {join, basename} = await import('path');
    const dapiDir = join(__dirname, 'dapi');
    const files = await readdir(dapiDir);
    for (const file of files) {
        if (file.endsWith('.js') && !testsExclude.includes(basename(file)))
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
    process.exit(0);
}


main().catch((err) => {
    console.error("Unhandled error:", err);
    process.exit(1);
});