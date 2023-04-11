import fs from 'fs';
import path from 'path';
import { exec } from 'child_process';
import * as utils from '../utils/utils';
import * as color from '../utils/color-utils';


const curDir = process.cwd();
const repositoryDir = path.dirname(path.dirname(curDir));

const apiPackageName = 'datagrok-api';
const libScope = '@datagrok-libraries';
const paths = {
  [apiPackageName]: path.join(repositoryDir, 'js-api'),
  [libScope]: path.join(repositoryDir, 'libraries'),
};

const packageDependencies: {[packagePath: string]: utils.Indexable} = {};

/** Links local packages. */
export function link(args: LinkArgs) {
  const nOptions = Object.keys(args).length - 1;
  if (nOptions > 1 || args['_'].length > 1 || (nOptions === 1 && !args.local && !args.npm))
    return false;
  const local = args.npm ? false : true;  // args.local is default

  if (!utils.isPackageDir(curDir)) {
    color.error('File `package.json` not found. Run the command from the package directory');
    return false;
  }

  for (const p of Object.values(paths)) {
    if (!fs.existsSync(p)) {
      color.error(`Directory ${p} not found. Run the command from the public package repository`);
      return false;
    }
  }

  const dependencies = readDependencies(curDir);

  if (local) {
    const packageHierarchy: string[] = getHierarchy(curDir);
    function link(packagePath: string): void {
      const packageJsonPath = path.join(packagePath, 'package.json');
      const json = JSON.parse(fs.readFileSync(packageJsonPath, 'utf-8'));
      json.scripts['link-all'] = generateLinkScript(packagePath, packageHierarchy);
      json.scripts['build-all'] = generateBuildScript(packagePath, packageHierarchy);
      fs.writeFileSync(packageJsonPath, JSON.stringify(json, null, 2), 'utf-8');
      exec('npm install && npm link && npm run link-all', { cwd: packagePath }, (err, stdout, stderr) => {
        if (err) throw err;
        else console.log(stderr, stdout);
      });
    }

    packageHierarchy.forEach(packageName => link(getPackagePath(packageName)));
    link(curDir);
  } else {
    runScript(curDir, 'npm install', dependencies, {
      dirMessage: 'Unlinking local packages in ',
      successMessage: 'Local packages have been successfully unlinked.',
    });
  }

  return true;
}

function readDependencies(packagePath: string): utils.Indexable {
  if (packagePath in packageDependencies)
    return packageDependencies[packagePath];
  const fileContent = fs.readFileSync(path.join(packagePath, 'package.json'), 'utf-8');
  const json = JSON.parse(fileContent);
  const libs: utils.Indexable = {};
  for (let dep in json.dependencies) {
    if (dep === apiPackageName || dep.startsWith(`${libScope}/`))
      libs[dep] = json.dependencies[dep];
  }
  packageDependencies[packagePath] = libs;
  return libs;
}

function getPackagePath(packageName: string): string {
  return packageName === apiPackageName ? paths[packageName] : path.join(paths[libScope], packageName.split('/')[1]);
}

/** Forms a hierarchy to understand in which order packages should be linked. */
function getHierarchy(packageDir: string): string[] {
  const hierarchy: string[] = [];
  const dependencies = Object.keys(readDependencies(packageDir));
  const cachedHierarchy: {[path: string]: string[]} = {};
  for (const dep of dependencies) {
    let idx = hierarchy.indexOf(dep);
    if (idx === -1)
      idx = hierarchy.push(dep) - 1;
    const depPath = getPackagePath(dep);
    const internalHierarchy = cachedHierarchy[depPath] ?? getHierarchy(depPath);
    if (!(depPath in cachedHierarchy))
      cachedHierarchy[depPath] = internalHierarchy;
    for (const internalDep of internalHierarchy) {
      const depIdx = hierarchy.indexOf(internalDep);
      if (depIdx === -1) {
        // Insert the internal dependency before the main package
        hierarchy.splice(idx, 0, internalDep);
        idx++;
      } else if (depIdx > idx) {
        // Remove the internal dependency from the list and place it before the main package
        // (doesn't affect the order of same-level libraries, their order is arbitrary)
        hierarchy.splice(depIdx, 1);
        idx = hierarchy.indexOf(dep);
        hierarchy.splice(idx, 0, internalDep);
        idx++;
      }
    }
  }
  return hierarchy;
}

/** Executes a script for a package and its dependencies with messages on the current progress. */
function runScript(packageDir: string, script: string, dependencies: utils.Indexable, options:
  {dirMessage?: string, successMessage?: string, callback?: (dir: string) => void} = {}) {
  exec(script, {cwd: packageDir}, (err, stdout, stderr) => {
    if (options.dirMessage)
      console.log(`${options.dirMessage}${path.basename(packageDir)}`);
    if (options.callback)
      options.callback(packageDir);
    if (err) throw err;
    else {
      console.log(stderr, stdout);
      for(const dep in dependencies) {
        const depPath = getPackagePath(dep);
        exec(script, {cwd: depPath}, (err, stdout, stderr) => {
          if (options.dirMessage)
            console.log(`${options.dirMessage}${path.basename(depPath)}`);
          if (options.callback)
            options.callback(depPath);
          if (err) throw err;
          else console.log(stderr, stdout);
        });
      }
      // if (options.successMessage)
      //   setTimeout(() => color.success(options.successMessage!), 5000);
    }
  });
}

/** Generates a package script to build all dependencies using the provided hierarchy. */
function generateBuildScript(packagePath: string, hierarchy: string[]): string {
  const dependencies = Object.keys(readDependencies(packagePath));
  const packageNames = hierarchy.filter((p) => dependencies.includes(p));
  const prefix = `./${path.relative(packagePath, repositoryDir).split(path.sep).join('/')}/`;
  let script = '';
  for (const packageName of packageNames) {
    script += `npm --prefix ${prefix}${packageName === apiPackageName ?
      'js-api' : `libraries/${packageName.split('/')[1]}`} run build && `;
  }
  return `${script ? script : ''}npm run build`;
}

/** Generates a package script to link all dependencies using the provided hierarchy. */
function generateLinkScript(packagePath: string, hierarchy: string[]): string {
  const dependencies = Object.keys(readDependencies(packagePath));
  const packageNames = hierarchy.filter((p) => dependencies.includes(p));
  for (const dep of dependencies) {
    if (!packageNames.includes(dep))
      color.error(`Hierarchy does not include package ${dep}`);
  }
  let script = `npm link${packageNames.length ? (' ' + packageNames.join(' ')) : ''}`;
  return script;
}

interface LinkArgs {
  _: string[],
  local?: boolean,
  npm?: boolean,
}
