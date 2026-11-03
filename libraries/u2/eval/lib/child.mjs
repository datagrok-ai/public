// node --import ./register.mjs child.mjs <grader.mjs> <solution.ts>
import {pathToFileURL} from 'node:url';
import {grade} from './harness.mjs';

const [graderPath, solutionPath] = process.argv.slice(2);
const {checks} = await import(pathToFileURL(graderPath).href);
process.stdout.write(`\n@@CHECKS ${JSON.stringify(Object.keys(checks))}\n`);
let mod;
let loadError;
try {
  mod = await import(pathToFileURL(solutionPath).href);
}
catch (e) {
  loadError = e ?? new Error('import failed');
}
await grade(checks, mod, loadError, {file: solutionPath});
process.exit(0);
