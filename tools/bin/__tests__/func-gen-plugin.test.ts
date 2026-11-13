import {describe, it, expect} from 'vitest';
import fs from 'fs';
import os from 'os';
import path from 'path';
import {fileURLToPath} from 'url';
import {spawnSync} from 'child_process';

const pluginPath = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..', '..', 'plugins', 'func-gen-plugin.js');

function makePackage(): string {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'grok-func-gen-'));
  fs.mkdirSync(path.join(dir, 'src'));
  fs.writeFileSync(path.join(dir, 'package.json'), JSON.stringify({name: '@datagrok/func-gen-test', friendlyName: 'FuncGenTest'}));
  fs.writeFileSync(path.join(dir, 'src', 'package.ts'), [
    'import * as grok from \'datagrok-api/grok\';',
    'export * from \'./package.g\';',
    'export class PackageFunctions {',
    '  @grok.decorators.func()',
    '  static async addedFunc(a: number, b: number): Promise<number> {',
    '    return a + b;',
    '  }',
    '}',
    '',
  ].join('\n'));
  fs.writeFileSync(path.join(dir, 'src', 'package.g.ts'), 'import * as DG from \'datagrok-api/dg\';\n');
  return dir;
}

/** Applies the plugin to a stub compiler and runs one compilation, in [dir] since `grok api` reads the cwd. */
function runPlugin(dir: string): void {
  const script = `
    const FuncGeneratorPlugin = require(${JSON.stringify(pluginPath)});
    const taps = [];
    const compiler = {context: process.cwd(), hooks: {compilation: {tap: (_name, fn) => taps.push(fn)}}};
    new FuncGeneratorPlugin({outputPath: './src/package.g.ts'}).apply(compiler);
    taps.forEach((fn) => fn({}));
  `;
  const res = spawnSync(process.execPath, ['-e', script], {cwd: dir, encoding: 'utf8'});
  expect(res.status, res.stderr).toBe(0);
}

describe('FuncGeneratorPlugin', () => {
  it('generates package-api.ts from the package.g.ts of the same build', () => {
    const dir = makePackage();
    try {
      runPlugin(dir);
      expect(fs.readFileSync(path.join(dir, 'src', 'package.g.ts'), 'utf8')).toContain('addedFunc');
      expect(fs.readFileSync(path.join(dir, 'src', 'package-api.ts'), 'utf8'))
        .toContain('grok.functions.call(\'FuncGenTest:AddedFunc\', { a, b })');
    } finally {
      fs.rmSync(dir, {recursive: true, force: true});
    }
  });
});
