import assert from 'node:assert/strict';
import {mkdirSync, mkdtempSync, rmSync, writeFileSync} from 'node:fs';
import {tmpdir} from 'node:os';
import {join, resolve} from 'node:path';
import {after, before, test} from 'node:test';
import {BASE_TIERS, findRoot, librarySource, loadProject, PACKAGE_NAME} from '../src/project.js';

let tmp: string;

before(() => {
  tmp = mkdtempSync(join(tmpdir(), 'bdd-project-'));
});

after(() => {
  rmSync(tmp, {recursive: true, force: true});
});

function pkg(name: string, config?: object, withBindings = true): string {
  const dir = join(tmp, name);
  mkdirSync(join(dir, 'bdd', 'features'), {recursive: true});
  if (withBindings)
    mkdirSync(join(dir, 'bdd', 'bindings'), {recursive: true});
  writeFileSync(join(dir, 'package.json'), JSON.stringify({name: `@datagrok/${name}`}));
  if (config)
    writeFileSync(join(dir, 'bdd', 'bdd.config.json'), JSON.stringify(config));
  return dir;
}

test('a package project lives under <package>/bdd', () => {
  const dir = pkg('demo');
  const project = loadProject(dir);
  assert.equal(resolve(project.root), resolve(join(dir, 'bdd')));
  assert.equal(resolve(project.packageDir), resolve(dir));
  assert.equal(project.name, '@datagrok/demo');
  assert.deepEqual(project.tiers, []);
  assert.deepEqual(project.sources.map((s) => s.dir.split(/[\\/]/).slice(-2).join('/')), ['bindings/common', 'bindings/platform', 'bdd/bindings']);
  assert.equal(resolve(project.featuresDir), resolve(join(dir, 'bdd', 'features')));
  assert.equal(resolve(project.generatedDir), resolve(join(dir, 'bdd', 'generated')));
});

test('a directory with features/ of its own is a project too', () => {
  const dir = join(tmp, 'flat');
  mkdirSync(join(dir, 'features'), {recursive: true});
  assert.equal(resolve(findRoot(dir)), resolve(dir));
  assert.equal(loadProject(dir).sources.length, BASE_TIERS.length);
});

test('tiers come from bdd.config.json and must exist in the library', () => {
  assert.throws(() => loadProject(pkg('bad-tier', {tiers: ['no-such-tier']})), /unknown tier "no-such-tier"/);
  assert.throws(() => loadProject(pkg('base-tier', {tiers: ['common']})), /always loaded/);
  const project = loadProject(pkg('viewers', {tiers: ['viewers']}, false));
  assert.deepEqual(project.tiers, ['viewers']);
  assert.ok(project.sources.some((s) => /tiers[\\/]viewers$/.test(s.dir)));
});

test('library modules are imported by package subpath, project modules by file', () => {
  const common = librarySource('common');
  assert.equal(common.specifierFor(join(common.dir, 'steps.ts')), `${PACKAGE_NAME}/bindings/common/steps`);
  const viewers = librarySource(join('tiers', 'viewers'));
  assert.equal(viewers.specifierFor(join(viewers.dir, 'steps.js')), `${PACKAGE_NAME}/bindings/tiers/viewers/steps`);
  const dir = pkg('local');
  const local = loadProject(dir).sources.at(-1)!;
  assert.equal(local.specifierFor(join(local.dir, 'elements.ts')), join(local.dir, 'elements.ts'));
});

test('no features, no project', () => {
  assert.throws(() => findRoot(join(tmp, 'nowhere')), /no bdd project/);
});
