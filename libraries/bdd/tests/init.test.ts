import assert from 'node:assert/strict';
import {existsSync, mkdirSync, mkdtempSync, readFileSync, rmSync, writeFileSync} from 'node:fs';
import {tmpdir} from 'node:os';
import {join} from 'node:path';
import {after, before, test} from 'node:test';
import {scaffold} from '../src/init.js';
import {PACKAGE_NAME} from '../src/project.js';

let tmp: string;

before(() => {
  tmp = mkdtempSync(join(tmpdir(), 'bdd-init-'));
});

after(() => {
  rmSync(tmp, {recursive: true, force: true});
});

function pkg(name: string, manifest: object): string {
  const dir = join(tmp, name);
  mkdirSync(dir, {recursive: true});
  writeFileSync(join(dir, 'package.json'), JSON.stringify(manifest, null, 2) + '\n');
  return dir;
}

test('init scaffolds bdd/, the editor settings, the ignore list and the manifest entries', () => {
  const dir = pkg('demo', {name: '@datagrok/demo', friendlyName: 'Demo App', scripts: {build: 'webpack'}});
  const result = scaffold(dir);
  for (const file of ['bdd/package.json', 'bdd/bdd.config.json', 'bdd/tsconfig.json', 'bdd/bindings/elements.ts',
    'bdd/bindings/steps.ts', 'bdd/features/smoke.feature', '.vscode/settings.json'])
    assert.ok(existsSync(join(dir, file)), `${file} exists`);
  assert.deepEqual(JSON.parse(readFileSync(join(dir, 'bdd', 'package.json'), 'utf8')), {type: 'module', private: true});
  assert.deepEqual(JSON.parse(readFileSync(join(dir, 'bdd', 'bdd.config.json'), 'utf8')), {tiers: []});
  assert.match(readFileSync(join(dir, 'bdd', 'bindings', 'elements.ts'), 'utf8'), /context\('Demo App app'/);
  assert.match(readFileSync(join(dir, 'bdd', 'bindings', 'steps.ts'), 'utf8'), /user opens the Demo App app/);
  assert.match(readFileSync(join(dir, 'bdd', 'features', 'smoke.feature'), 'utf8'), /Given user is logged in/);
  assert.match(readFileSync(join(dir, '.gitignore'), 'utf8'), /bdd\/test-results\//);
  const manifest = JSON.parse(readFileSync(join(dir, 'package.json'), 'utf8'));
  assert.ok(manifest.devDependencies[PACKAGE_NAME]);
  assert.ok(manifest.devDependencies['@playwright/test']);
  assert.equal(manifest.scripts['test:bdd'], 'grok-bdd run');
  assert.equal(manifest.scripts.build, 'webpack');
  assert.ok(result.created.includes('bdd/features/smoke.feature'));
  assert.equal(result.skipped.length, 0);
});

test('init never overwrites: a second run reports everything as existing', () => {
  const dir = pkg('again', {name: 'again'});
  scaffold(dir);
  writeFileSync(join(dir, 'bdd', 'features', 'smoke.feature'), 'Feature: mine\n');
  const before = readFileSync(join(dir, 'package.json'), 'utf8');
  const result = scaffold(dir);
  assert.deepEqual(result.created, []);
  assert.ok(result.skipped.some((s) => s.startsWith('bdd/features')));
  assert.ok(result.skipped.includes('package.json'));
  assert.equal(readFileSync(join(dir, 'bdd', 'features', 'smoke.feature'), 'utf8'), 'Feature: mine\n');
  assert.equal(readFileSync(join(dir, 'package.json'), 'utf8'), before);
});

test('init merges missing cucumber keys into an existing settings file and keeps the rest', () => {
  const dir = pkg('settings', {name: 'settings'});
  mkdirSync(join(dir, '.vscode'));
  writeFileSync(join(dir, '.vscode', 'settings.json'), JSON.stringify({'editor.tabSize': 2}));
  scaffold(dir);
  const settings = JSON.parse(readFileSync(join(dir, '.vscode', 'settings.json'), 'utf8'));
  assert.equal(settings['editor.tabSize'], 2);
  assert.deepEqual(settings['cucumber.features'], ['bdd/features/**/*.feature']);
});

test('init adds no smoke feature to a project that has features of its own', () => {
  const dir = pkg('own', {name: 'own'});
  mkdirSync(join(dir, 'bdd', 'features', 'deep'), {recursive: true});
  writeFileSync(join(dir, 'bdd', 'features', 'deep', 'mine.feature'), 'Feature: mine\n');
  const result = scaffold(dir);
  assert.ok(!existsSync(join(dir, 'bdd', 'features', 'smoke.feature')));
  assert.ok(result.skipped.some((s) => s.startsWith('bdd/features')));
});

test('init refuses a directory without package.json', () => {
  const dir = join(tmp, 'nowhere');
  mkdirSync(dir);
  assert.throws(() => scaffold(dir), /not a package/);
});
