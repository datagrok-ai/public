
/* eslint-disable valid-jsdoc */
import fs from 'fs';
import path from 'path';
import { runScript } from '../utils/utils';

const repositoryDirNameRegex = new RegExp(path.join('1', '2')[1] + 'public$');

const excludedPackages: string[] = ['@datagrok/diff-grok', ''];

const curDir = process.cwd();
let repositoryDir = curDir;
let currentPackage: PackageData | undefined;
let localPackageDependencies: PackageData[];
let packagesToLink: Set<string>;
while (path.dirname(repositoryDir) !== repositoryDir) {
  if (repositoryDirNameRegex.test(repositoryDir))
    break;
  repositoryDir = path.dirname(repositoryDir);
}

let verbose = false;
let pathMode = false;
let devMode = false;
let unlink = false;

const apiDir = path.join(repositoryDir, 'js-api');
const libDir = path.join(repositoryDir, 'libraries');
const packageDir = path.join(repositoryDir, 'packages');

const apiPackageName = 'datagrok-api';
const libName = '@datagrok-libraries/';
const packageName = '@datagrok/';

export async function link(args: LinkArgs) {
  verbose = args.verbose ?? false;
  devMode = args.dev ?? false;
  pathMode = args.path ?? false;
  unlink = args.unlink ?? false;

  localPackageDependencies = []
  packagesToLink = new Set<string>();
  collectPackagesData(curDir);
  currentPackage = new PackageData(curDir);

  if (unlink) {
    await unlinkPackages();
    console.log('Package unlinked')
  }
  else {
    await linkPackages();
    if (pathMode)
      console.log('Updated dependencies to local in package.json')
    else
      console.log('All packages/libraries linked to your package by nmp link command')
  }

  return true;
}

function collectPackagesData(packagePath: string = curDir): { dependencies: string[], version: string, name: string } {
  const packageJsonPath = path.join(packagePath, 'package.json');
  if (!fs.existsSync(packageJsonPath))
    return { dependencies: [], version: '0.0.1', name: '' };
  const json = JSON.parse(fs.readFileSync(packageJsonPath, 'utf-8'));
  let dependencies: string[] = [];
  dependencies = dependencies.concat(collectPacakgeDataFromJsonObject(json.dependencies));
  if (devMode)
    dependencies = dependencies.concat(collectPacakgeDataFromJsonObject(json.devDependencies));
  return { dependencies: dependencies, version: json.version, name: json.name };
}

function collectPacakgeDataFromJsonObject(object: any): string[] {
  let result: string[] = [];
  for (let dependencyName of Object.keys(object)) {
    let nameSplitted = dependencyName.split('/');
    if (dependencyName === apiPackageName)
      result = result.concat(parsePackageDependencies(dependencyName, apiDir));
    else if (dependencyName.includes(libName))
      result = result.concat(parsePackageDependencies(dependencyName, path.join(libDir, nameSplitted[nameSplitted.length - 1])));
    else if (dependencyName.includes(packageName))
      result = result.concat(parsePackageDependencies(dependencyName, path.join(packageDir, toCamelCase(nameSplitted[nameSplitted.length - 1]))));
  }
  return result;
}

function toCamelCase(input: string): string {
  return input
    .split('-')
    .map(word => word.charAt(0).toUpperCase() + word.slice(1))
    .join('');
}

function parsePackageDependencies(dependencyName: string, pathToLink: string): string[] {
  let result: string[] = [];
  if (!packagesToLink.has(dependencyName)) {
    packagesToLink.add(dependencyName);
    localPackageDependencies.push(new PackageData(pathToLink));
  }
  result.push(dependencyName);
  return result;
}

async function unlinkPackages() {
  const packaegs = [...localPackageDependencies];
  if (currentPackage)
    packaegs.push(currentPackage);
  for (let packageData of packaegs) {
    if (excludedPackages.includes(packageData.name))
      continue;

    if (verbose)
      console.log(`Package ${packageData.name} is unlinking`);

    const packageJsonPath = path.join(packageData.packagePath, 'package.json');
    const packageJson = JSON.parse(fs.readFileSync(packageJsonPath, 'utf8'));
    packageJson.dependencies = updateDependenciesToVersion(packageData, packageJson.dependencies);
    if (devMode)
      packageJson.devDependencies = updateDependenciesToVersion(packageData, packageJson.devDependencies);
    fs.writeFileSync(packageJsonPath, JSON.stringify(packageJson, null, 2));

    await runScript(`npm unlink ${packageData.name}`, packageData.packagePath);
  }
}

function updateDependenciesToVersion(packageData: PackageData, dependecyNode: any) {
  for (let dependency of packageData.dependencies) {
    if (excludedPackages.includes(dependency))
      continue;
    if (dependecyNode[dependency]) {
      const packageToLink = (localPackageDependencies.filter((e) => e.name === dependency) ?? [])[0];
      if (packageToLink) {
        const startIndex = packageToLink.packagePath.indexOf(`public${path.sep}`) + `public${path.sep}`.length;
        if (startIndex !== -1)
          dependecyNode[dependency] = `^${packageToLink.version}`;
      }
    }
  }
  return dependecyNode;
}

async function linkPackages() {
  let anyChanges = true;
  for (let element of excludedPackages)
    packagesToLink.delete(element);

  if (verbose) {
    console.log('Packages to link:')
    console.log(localPackageDependencies);
  }

  while (anyChanges && packagesToLink.size > 0) {
    anyChanges = false;
    let mapElements = localPackageDependencies.filter(x => x.dependencies.every(i => !packagesToLink.has(i)) && packagesToLink.has(x.name) && !excludedPackages.includes(x.name));
    if (mapElements.length === 0)
      break;
    for (let element of mapElements) {

      if (verbose)
        console.log(`Package ${element.name} is linking`)
      if (pathMode)
        await linkPathMode(element);
      await runScript(`npm install`, element.packagePath);
      if (pathMode)
        await runScript(`npm run build`, element.packagePath);

      if (!pathMode)
        await linkNpmMode(element);
      packagesToLink.delete(element.name);
      anyChanges = true;
    }

    if (anyChanges === false)
      throw (new Error(`There is loop with next packages: ${JSON.stringify(Array.from(packagesToLink)).toString()}`));
  }

  let names = localPackageDependencies.map(x => x.name);

  if (currentPackage && pathMode)
    await linkPathMode(currentPackage);
  await runScript(`npm install`, curDir);
  if (pathMode === false)
    await runScript(`npm link ${names.join(' ')}`, curDir);
}

async function linkNpmMode(packageData: PackageData) {
  await runScript(`npm link ${packageData.dependencies.join(' ')}`, packageData.packagePath);
  await runScript(`npm link`, packageData.packagePath);
}

async function linkPathMode(packageData: PackageData) {
  const packageJsonPath = path.join(packageData.packagePath, 'package.json');
  const packageJson = JSON.parse(fs.readFileSync(packageJsonPath, 'utf8'));
  packageJson.dependencies = packageJson.dependencies || {};
  packageJson.devDependencies = packageJson.devDependencies || {};
  packageJson.dependencies = updateDependenciesToLocal(packageData, packageJson.dependencies);
  if (devMode)
    packageJson.devDependencies = updateDependenciesToLocal(packageData, packageJson.devDependencies);

  fs.writeFileSync(packageJsonPath, JSON.stringify(packageJson, null, 2));
}

function updateDependenciesToLocal(packageData: PackageData, dependecyNode: any) {
  const dirNames = packageData.packagePath.split(path.sep);
  const backRoute = `../`.repeat(dirNames.length - dirNames.indexOf('public') - 1);

  for (let dependency of packageData.dependencies) {
    if (excludedPackages.includes(dependency))
      continue;
    if (dependecyNode[dependency]) {
      const packageToLink = (localPackageDependencies.filter((e) => e.name === dependency) ?? [])[0];
      if (packageToLink) {
        const startIndex = packageToLink.packagePath.indexOf(`public${path.sep}`) + `public${path.sep}`.length;
        if (startIndex !== -1)
          dependecyNode[dependency] = `${backRoute}${packageToLink.packagePath.substring(startIndex).replaceAll(path.sep, '/')}`;
      }
    }
  }
  return dependecyNode;
}

class PackageData {
  name: string;
  packagePath: string;
  dependencies: string[] = [];
  version: string;

  constructor(packagePath: string) {
    let packageJsonData = collectPackagesData(packagePath);
    this.name = packageJsonData.name;
    this.packagePath = packagePath;
    this.dependencies = packageJsonData.dependencies;
    this.version = packageJsonData.version;
  }
}

interface LinkArgs {
  verbose?: boolean,
  unlink?: boolean,
  path?: boolean,
  dev?: boolean,
}