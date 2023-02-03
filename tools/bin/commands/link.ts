import fs from 'fs';
import path from 'path';
import { exec } from 'child_process';
import * as utils from '../utils/utils';
import * as color from '../utils/color-utils';

/**
 * Link: api > utils > other libs - npm install, npm run build, npm link
 * Unlink: npm install in the package and its dependencies, npm unlink
 */

const apiPackageName = 'datagrok-api';
const libScope = '@datagrok-libraries';
const paths = {
  [apiPackageName]: path.join(path.dirname(path.dirname(path.dirname(__dirname))), 'js-api'),
  [libScope]: path.join(path.dirname(path.dirname(path.dirname(__dirname))), 'libraries'),
};

/** Links local packages. */
export function link(args: LinkArgs) {
  const nOptions = Object.keys(args).length - 1;
  if (nOptions > 0 || args['_'].length > 1)
    return false;
  const curDir = process.cwd();
  const packageDir = path.join(curDir, 'package.json');

  if (!utils.isPackageDir(curDir)) {
    color.error('File `package.json` not found. Run the command from the package directory');
    return false;
  }

  const dependencies = readDependencies(JSON.parse(fs.readFileSync(packageDir, 'utf-8')));
  const modulesDir = path.join(curDir, 'node_modules');
  if (!fs.existsSync(modulesDir)) {
    console.log('Running `npm install` to get the required dependencies...\n');
    exec('npm install', (err, stdout, stderr) => {
      if (err) throw err;
      else console.log(stderr, stdout);
    });
  }
  // The order should start with js-api, then libraries/utils, then other libraries
  let apiModule = null;
  let hasUtils = false;
  let libs: fs.Dirent[] = [];
  const allModules = fs.readdirSync(modulesDir, {encoding: 'utf-8', withFileTypes: true});
  for (const m of allModules) {
    if (m.name === apiPackageName)
      apiModule = m;
    else if (m.name === libScope) {
      libs = fs.readdirSync(path.join(modulesDir, m.name), {encoding: 'utf-8', withFileTypes: true});
      for (const lib of libs) {
        if (lib.name === 'utils')
          hasUtils = true;
      }
    }
  }

  function installLibs(libs: fs.Dirent[]) {
    if (!libs.length)
      return;
    const f = () => {
      for (const lib of libs) {
        if (`${libScope}/${lib.name}` in dependencies) {
          if (lib.isSymbolicLink()) {
            console.log(`Package "${lib.name}" is linked. Skipping...`);
          } else {
            console.log(`Linking "${lib.name}"...`);
            exec('npm install && npm run build && npm link',  { cwd: path.join(paths[libScope], lib.name) }, (err, stdout, stderr) => {
              if (err) throw err;
              else console.log(stderr, stdout);
            });
          }
        }
      }
    };
    if (hasUtils) {
      const utilsModule = libs.splice(libs.findIndex((m) => m.name === 'utils'), 1)[0];
      if (utilsModule.isSymbolicLink()) {
        console.log(`Package "${utilsModule.name}" is linked. Skipping...`);
      } else {
        console.log(`Linking "${utilsModule.name}"...`);
        exec('npm install && npm run build && npm link',  { cwd: path.join(paths[libScope], 'utils') }, (err, stdout, stderr) => {
          if (err) throw err;
          else {
            console.log(stderr, stdout);
            f();
          }
        });
      }
    } else {
      f();
    }
  }

  if (apiModule != null) {
    if (apiModule.isSymbolicLink()) {
      console.log(`Package "${apiModule.name}" is linked. Skipping...`);
      installLibs(libs);
    } else {
      console.log(`Linking "${apiModule.name}"...`);
      exec(`npm install && npm run build && npm link`, { cwd: paths[apiPackageName] }, (err, stdout, stderr) => {
        if (err) throw err;
        else {
          console.log(stderr, stdout);
          installLibs(libs);
          const packageNames = Object.keys(dependencies).join(' ');
          exec(`npm link ${packageNames}`);
        }
      });
    }
  }

  return true;
}

/** Unlinks local packages and runs `npm i`. */
export function unlink(args: LinkArgs) {
  const nOptions = Object.keys(args).length - 1;
  if (nOptions > 0 || args['_'].length > 1)
    return false;

  const curDir = process.cwd();
  const packageDir = path.join(curDir, 'package.json');

  if (!utils.isPackageDir(curDir)) {
    color.error('File `package.json` not found. Run the command from the package directory');
    return false;
  }

  const dependencies = readDependencies(JSON.parse(fs.readFileSync(packageDir, 'utf-8')));
  const packageNames = Object.keys(dependencies).join(' '); // npm unlink ${packageNames}
  exec(`npm install`, (err, stdout, stderr) => {
    console.log(`Unlinking local packages in ${path.basename(curDir)}`);
    if (err) throw err;
    else {
      console.log(stderr, stdout);
      for (const dep in dependencies) {
        const depPath = dep === apiPackageName ? paths[dep] : path.join(paths[libScope], dep.split('/')[1]);
        exec(`npm install`, {cwd: depPath}, (err, stdout, stderr) => {
          console.log(`Unlinking local packages in ${path.basename(depPath)}`);
          if (err) throw err;
          else console.log(stderr, stdout);
        });
      }
    }
  });

  return true;
}

function readDependencies(json: utils.Indexable): utils.Indexable {
  const libs: utils.Indexable = {};
  for (let dep in json.dependencies) {
    if (dep === apiPackageName || dep.startsWith(`${libScope}/`))
      libs[dep] = json.dependencies[dep];
  }
  return libs;
}

interface LinkArgs {
  _: string[],
}
