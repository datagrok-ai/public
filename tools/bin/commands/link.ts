
/* eslint-disable valid-jsdoc */
import fs from 'fs';
import path from 'path';
import { runScript } from '../utils/utils';
import { glob } from 'glob';

const excludedPackages: string[] = ['@datagrok-misc/eslint-plugin-config', 'wiki-merge', 'datagrok-tools', ''];
const versionDependencyRegex = /(?:\^)?\d+\.\d+\.\d+/;
let packagesInRepo: Record<string, { path: string, version: string }>;
let packagesOutOfRepo: Record<string, { path: string, version: string }>;


const curDir = process.cwd();
let repositoryDir = curDir;
let containerDir = curDir;
let dirStep = path.dirname(curDir);
let currentPackage: PackageData | undefined;
let localPackageDependencies: PackageData[];
let packagesToLink: Set<string>;

while (path.dirname(dirStep) !== dirStep) {
  if (fs.existsSync(path.join(dirStep, '.git')) && repositoryDir === curDir)
    repositoryDir = dirStep;
  if (!fs.existsSync(path.join(dirStep, '.git')) && repositoryDir !== curDir && containerDir === curDir) {
    containerDir = dirStep;
    break;
  }
  dirStep = path.dirname(dirStep);
}

let verbose = false;
let pathMode = false;
let devMode = false;
let unlink = false;

export async function link(args: LinkArgs) {
  verbose = args.verbose ?? false;
  devMode = args.dev ?? false;
  pathMode = args.path ?? false;
  unlink = args.unlink ?? false; 

  let collectedPackages = collectAvaliablePackages(args['without-common-dir']);
  packagesInRepo = collectedPackages.packagesInRepo;
  packagesOutOfRepo = collectedPackages.packagesOutOfRepo;

  localPackageDependencies = []
  packagesToLink = new Set<string>();
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
  for (let dependencyName of Object.keys(object ?? {})) {
    if (packagesInRepo[dependencyName])
      result = result.concat(parsePackageDependencies(dependencyName, path.dirname(packagesInRepo[dependencyName].path)));
    else if (packagesOutOfRepo[dependencyName])
      result = result.concat(parsePackageDependencies(dependencyName, path.dirname(packagesOutOfRepo[dependencyName].path)));

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

function collectAvaliablePackages(noOutLink: boolean = false): { packagesInRepo: Record<string, { path: string, version: string }>, packagesOutOfRepo: Record<string, { path: string, version: string }> } {
  let repositoryPackages = collectAvaliablePackagesPathesFromDir(repositoryDir);
  let commonPakcages = collectAvaliablePackagesPathesFromDir(containerDir, repositoryDir);

  const packagesInRepoBaseInfo = parsePackages(repositoryPackages);
  const packagesOutOfRepoBaseInfo = noOutLink? [] : parsePackages(commonPakcages);
  const packagesInRepo = Object.fromEntries(packagesInRepoBaseInfo.map(({ name, ...rest }) => [name, rest]));
  const packagesOutOfRepo = Object.fromEntries(packagesOutOfRepoBaseInfo.map(({ name, ...rest }) => [name, rest]));
  return { packagesInRepo: packagesInRepo, packagesOutOfRepo: packagesOutOfRepo };
}

function collectAvaliablePackagesPathesFromDir(dir: string, dirToExclude: string = '') {
  let res = glob.sync('./**/package.json', {
    cwd: dir,
    ignore: ['**/node_modules/**'],
  }).map((e) => path.join(dir, e));
  if (dirToExclude.length > 0)
    res = res.filter((e) => !e.includes(dirToExclude))
  return res;
}

function parsePackages(packagPathes: string[]) {
  let res = packagPathes.map((e) => {
    let packageData: any = {};
    if (fs.existsSync(e)) {
      packageData.path = e;
      const packageJson = JSON.parse(fs.readFileSync(e, 'utf8'));
      packageData.name = packageJson.name;
      packageData.version = packageJson.version;
    }
    return packageData;
  })
  return res;
}

async function unlinkPackages() {
  const packages = [...localPackageDependencies];
  if (currentPackage)
    packages.push(currentPackage);
  for (let packageData of packages) {
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
    if (dependecyNode[dependency] && !versionDependencyRegex.test(dependecyNode[dependency])) {
      const packageToLink = (localPackageDependencies.filter((e) => e.name === dependency) ?? [])[0];
      if (packageToLink) {
        dependecyNode[dependency] = `^${packagesInRepo[dependency]?.version ?? packagesOutOfRepo[dependency]?.version}`;
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

  const backRouteForCont = `../`.repeat(packageData.packagePath.replace(containerDir, '').split(path.sep).length - 1);
  const backRouteForRepo = `../`.repeat(packageData.packagePath.replace(repositoryDir, '').split(path.sep).length - 1);

  for (let dependency of packageData.dependencies) {
    if (excludedPackages.includes(dependency))
      continue;
    if (dependecyNode[dependency]) {
      const packageToLink = (localPackageDependencies.filter((e) => e.name === dependency) ?? [])[0];
      if (packageToLink) {
        if (packageToLink.packagePath.includes(repositoryDir) && packageData.packagePath.includes(repositoryDir))
          dependecyNode[dependency] = `${backRouteForRepo}${packageToLink.packagePath.replace(repositoryDir, '').replaceAll(path.sep, '/').replace('/', '')}`;
        else
          dependecyNode[dependency] = `${backRouteForCont}${packageToLink.packagePath.replace(containerDir, '').replaceAll(path.sep, '/').replace('/', '')}`;

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
  'without-common-dir'?: boolean,
  path?: boolean,
  dev?: boolean,
}