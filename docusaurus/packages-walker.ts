import * as fs from 'fs';
import * as path from 'path';
import { exec } from 'child_process';
import { promisify } from 'util';

//invocation: 
//ts-node .\packages-walker.ts --commands="command1 && command2 ..." --rootDir="direcotry"
//without rootDir input script uses package dir of the public repo
//without commands input script runs npm i command

const execAsync = promisify(exec);

const argv = require('minimist')(process.argv.slice(2), {
  alias: { k: 'key', h: 'help', r: 'recursive' },
});

const rootDir = argv.rootDir ?? path.join(...[process.cwd(), '..', 'packages']);
const commands = argv.commands?.split('&&') ?? ['npm install'];
const verbose = argv.verbose;

if (!fs.existsSync(rootDir)) {
  console.error(`Root directory not found: ${rootDir}`);
  process.exit(1);
}

if (commands.length === 0) {
  console.error(`No commands provided to invoke`);
  process.exit(1);
}

async function delay(ms: number): Promise<void> {
  return new Promise(resolve => setTimeout(resolve, ms));
}

async function walkAllPackages(): Promise<void> {
  const packages = fs.readdirSync(rootDir);
  for (const packageName of packages) {
    const packagePath = path.join(...[rootDir, packageName]);
    const packageJsonPath = path.join(...[rootDir, packageName, 'package.json']);
    if (fs.statSync(packagePath).isDirectory() && fs.existsSync(packageJsonPath)) {
      const nodeModulesDir = path.join(...[rootDir, packageName, 'node_modules']);
      // if (fs.existsSync(nodeModulesDir))
      //   fs.rmdirSync(nodeModulesDir, {recursive:true});
      for (const command of commands) {
        console.log(`${packagePath}: ${command} - Started`);
        if (!(await runCommand(command, packagePath))) {
          console.log(`${packagePath}: ${command} - Error`);
          break;
        }
        else
          console.log(`${packagePath}: ${command} - Finished`);
        await delay(1000);
      }
    }
  }
  return;
}

async function runCommand(command: string, directory: string): Promise<boolean> {
  const { stdout, stderr } = await execAsync(command, { cwd: directory });
  let result = true;
  if (verbose)
    console.log(stdout)
  if (stderr?.length > 0 && stderr.includes('npm error')) {
    if (verbose)
      console.log(stderr)
    result = false;
  }
  return true;
}

(async () => {
  await walkAllPackages();
})();