/* `grok-bdd init`: the bdd/ project of a package — folders, config, the bindings to start from, a
   smoke feature that passes on any stand, the manifest entries, the editor settings. Run in the
   package directory (where package.json is). Never overwrites: what exists is reported as such. */
import {appendFileSync, existsSync, mkdirSync, readFileSync, writeFileSync} from 'node:fs';
import {join, sep} from 'node:path';
import {listFiles} from './discover.js';
import {LIB_ROOT, PACKAGE_NAME} from './project.js';

export interface InitResult {
  created: string[];
  skipped: string[];
  notes: string[];
}

interface Manifest {
  name?: string;
  friendlyName?: string;
  scripts?: Record<string, string>;
  devDependencies?: Record<string, string>;
  [key: string]: unknown;
}

const GITIGNORE = ['bdd/test-results/', 'bdd/e2e/', 'bdd/.auth.json'];
const STATES = 'visible|hidden|present|absent|enabled|disabled|checked|unchecked|selected|empty|expanded|collapsed|focused';

function libraryVersions(): {bdd: string; playwright: string} {
  const lib = JSON.parse(readFileSync(join(LIB_ROOT, 'package.json'), 'utf8')) as Manifest;
  return {bdd: `^${lib.version as string}`, playwright: lib.devDependencies?.['@playwright/test'] ?? '>=1.40.0'};
}

function indentOf(json: string): string {
  return /^( +|\t)"/m.exec(json)?.[1] ?? '  ';
}

function titleOf(manifest: Manifest, short: string): string {
  return manifest.friendlyName?.trim() || short.charAt(0).toUpperCase() + short.slice(1);
}

function elementsTemplate(title: string): string {
  return `/* The names this package's features use. Platform names (toolbox, browse tab, context panel,
   console, status bar, open tableview, grid) are reserved; the app's own names live on a context,
   and apply after a step declared with \`enters\` (see steps.ts). Most of a u2 page needs no entry:
   "run button in toolbar" or "name input in dialog" resolve from the u2 contract alone. */
import {context} from '${PACKAGE_NAME}';

export const app = context('${title} app', {selector: '[data-u2-name="${title}"]'});
// app.element('results', {selector: '[data-u2-name="results"]', aliases: ['results panel']});
`;
}

function stepsTemplate(title: string, short: string): string {
  return `/* The steps only this package can define: how to reach its app. Everything generic (clicks,
   typing, selects, assertions) is the library's — \`grok-bdd list-steps\` prints it all. */
import type {Page} from '@playwright/test';
import {Given} from '${PACKAGE_NAME}';

/** The app's route on the stand — adjust: one app is \`/apps/<Package>\`, several are
 * \`/apps/<Package>/<App>\`. */
export const APP_PATH = '/apps/${short}';

export const openApp = Given('user opens the ${title} app', async (page: Page) => {
  await page.goto(APP_PATH, {waitUntil: 'domcontentloaded'});
  await page.locator('[data-u2-name="${title}"]').waitFor({timeout: 60000});
}, {tier: 'ui', enters: '${title} app', description: 'needs this package published on the stand'});
`;
}

function featureTemplate(title: string): string {
  return `Feature: ${title} smoke
  The first feature of this package: platform steps only, so it passes on any stand. Replace it
  with the app's own scenarios — "user opens the ${title} app" (bindings/steps.ts) is the way in.

  Scenario: The platform is up
    Given user is logged in
    Then browse tab should be visible
`;
}

const TSCONFIG = `{
  "compilerOptions": {
    "target": "es2022",
    "module": "NodeNext",
    "moduleResolution": "NodeNext",
    "lib": ["es2022", "dom"],
    "strict": true,
    "noEmit": true,
    "skipLibCheck": true,
    "esModuleInterop": true
  },
  "include": ["bindings", "generated"]
}
`;

const VSCODE = {
  'cucumber.features': ['bdd/features/**/*.feature'],
  'cucumber.glue': ['bdd/bindings/**/*.ts', `node_modules/${PACKAGE_NAME}/bindings/**/*.ts`],
  'cucumber.parameterTypes': [{name: 'state', regexp: STATES}],
};

export function scaffold(packageDir: string): InitResult {
  const manifestFile = join(packageDir, 'package.json');
  if (!existsSync(manifestFile))
    throw new Error(`${packageDir} is not a package: no package.json here — run \`grok-bdd init\` in the package directory`);
  const manifestText = readFileSync(manifestFile, 'utf8');
  const manifest = JSON.parse(manifestText) as Manifest;
  const short = (manifest.name ?? 'app').replace(/^@[^/]+\//, '');
  const title = titleOf(manifest, short);
  const result: InitResult = {created: [], skipped: [], notes: []};
  const rel = (file: string) => file.slice(packageDir.length + 1).split(sep).join('/');
  const write = (file: string, content: string) => {
    if (existsSync(file)) {
      result.skipped.push(rel(file));
      return;
    }
    mkdirSync(join(file, '..'), {recursive: true});
    writeFileSync(file, content, 'utf8');
    result.created.push(rel(file));
  };

  const bdd = join(packageDir, 'bdd');
  write(join(bdd, 'package.json'), '{"type": "module", "private": true}\n');
  write(join(bdd, 'bdd.config.json'), '{"tiers": []}\n');
  write(join(bdd, 'tsconfig.json'), TSCONFIG);
  write(join(bdd, 'bindings', 'elements.ts'), elementsTemplate(title));
  write(join(bdd, 'bindings', 'steps.ts'), stepsTemplate(title, short));
  // the smoke feature is a first feature, not one more: a project with features of its own keeps them
  if (listFiles(join(bdd, 'features'), '.feature').length > 0)
    result.skipped.push('bdd/features (has features already)');
  else
    write(join(bdd, 'features', 'smoke.feature'), featureTemplate(title));

  const settingsFile = join(packageDir, '.vscode', 'settings.json');
  if (!existsSync(settingsFile))
    write(settingsFile, JSON.stringify(VSCODE, null, 2) + '\n');
  else {
    try {
      const settings = JSON.parse(readFileSync(settingsFile, 'utf8')) as Record<string, unknown>;
      const missing = Object.entries(VSCODE).filter(([key]) => !(key in settings));
      if (missing.length === 0)
        result.skipped.push(rel(settingsFile));
      else {
        for (const [key, value] of missing)
          settings[key] = value;
        writeFileSync(settingsFile, JSON.stringify(settings, null, 2) + '\n', 'utf8');
        result.created.push(`${rel(settingsFile)} (${missing.map(([key]) => key).join(', ')})`);
      }
    } catch {
      result.notes.push(`${rel(settingsFile)} is not plain JSON — add the cucumber.features/glue/parameterTypes settings by hand (see the README)`);
    }
  }

  const gitignore = join(packageDir, '.gitignore');
  const ignored = existsSync(gitignore) ? readFileSync(gitignore, 'utf8') : '';
  const lines = GITIGNORE.filter((line) => !ignored.split(/\r?\n/).includes(line));
  if (lines.length > 0) {
    const eol = ignored.length === 0 || ignored.endsWith('\n') ? '' : '\n';
    appendFileSync(gitignore, `${eol}${ignored.length === 0 ? '' : '\n'}# bdd (grok-bdd) run artifacts\n${lines.join('\n')}\n`);
    result.created.push(`.gitignore (${lines.join(', ')})`);
  }
  else
    result.skipped.push('.gitignore');

  const versions = libraryVersions();
  const changes: string[] = [];
  manifest.devDependencies ??= {};
  for (const [name, version] of [[PACKAGE_NAME, versions.bdd], ['@playwright/test', versions.playwright]]) {
    if (!manifest.devDependencies[name]) {
      manifest.devDependencies[name] = version;
      changes.push(`devDependencies.${name}`);
    }
  }
  manifest.scripts ??= {};
  if (!manifest.scripts['test:bdd']) {
    manifest.scripts['test:bdd'] = 'grok-bdd run';
    changes.push('scripts.test:bdd');
  }
  if (changes.length > 0) {
    writeFileSync(manifestFile, JSON.stringify(manifest, null, indentOf(manifestText)) + '\n', 'utf8');
    result.created.push(`package.json (${changes.join(', ')})`);
  }
  else
    result.skipped.push('package.json');
  return result;
}
