import * as fs from 'fs';
import * as os from 'os';
import * as path from 'path';
import {spawnSync} from 'child_process';

const BIN = path.resolve(__dirname, '..', 'dist', 'cli.js');

function runCli(args: string[], cwd: string): {stdout: string; stderr: string; code: number} {
  const r = spawnSync('node', [BIN, ...args], {cwd, encoding: 'utf8'});
  return {stdout: r.stdout ?? '', stderr: r.stderr ?? '', code: r.status ?? -1};
}

function mkTmp(): string {
  return fs.mkdtempSync(path.join(os.tmpdir(), 'codegen-cli-'));
}

function rmTmp(dir: string): void {
  fs.rmSync(dir, {recursive: true, force: true});
}

const SRC_FOO = `// @async-source: foo-sync.ts
export async function fooAsync(x: number): Promise<number> { return await Promise.resolve(x); }
`;

beforeAll(() => {
  if (!fs.existsSync(BIN))
    throw new Error(`CLI not built — run 'npm run build' first. Expected: ${BIN}`);
});

describe('CLI', () => {
  test('--roots <dir>: finds and writes sync sibling', () => {
    const tmp = mkTmp();
    try {
      fs.mkdirSync(path.join(tmp, 'src'));
      fs.writeFileSync(path.join(tmp, 'src', 'foo.ts'), SRC_FOO);
      const {code, stdout} = runCli(['--roots', 'src'], tmp);
      expect(code).toBe(0);
      expect(stdout).toMatch(/wrote .*foo-sync\.ts/);
      const generated = fs.readFileSync(path.join(tmp, 'src', 'foo-sync.ts'), 'utf8');
      expect(generated).toContain('function fooAsync(');
      expect(generated).not.toContain('await');
    } finally { rmTmp(tmp); }
  });

  test('--roots accepts multiple directories', () => {
    const tmp = mkTmp();
    try {
      fs.mkdirSync(path.join(tmp, 'a'));
      fs.mkdirSync(path.join(tmp, 'b'));
      fs.writeFileSync(path.join(tmp, 'a', 'foo.ts'), SRC_FOO);
      fs.writeFileSync(path.join(tmp, 'b', 'bar.ts'),
        `// @async-source: bar-sync.ts\nexport async function barAsync(): Promise<number> { return 2; }\n`);
      const {code} = runCli(['--roots', 'a', 'b'], tmp);
      expect(code).toBe(0);
      expect(fs.existsSync(path.join(tmp, 'a', 'foo-sync.ts'))).toBe(true);
      expect(fs.existsSync(path.join(tmp, 'b', 'bar-sync.ts'))).toBe(true);
    } finally { rmTmp(tmp); }
  });

  test('--check: exits 0 when sync mirror is up-to-date', () => {
    const tmp = mkTmp();
    try {
      fs.mkdirSync(path.join(tmp, 'src'));
      fs.writeFileSync(path.join(tmp, 'src', 'foo.ts'), SRC_FOO);
      // First seed the sync sibling
      const seed = runCli(['--roots', 'src'], tmp);
      expect(seed.code).toBe(0);
      // Now check should be clean
      const {code, stdout, stderr} = runCli(['--check', '--roots', 'src'], tmp);
      expect(code).toBe(0);
      expect(stderr).toBe('');
      expect(stdout).toMatch(/up to date/);
    } finally { rmTmp(tmp); }
  });

  test('--check: exits 1 when sync mirror drifted; does not rewrite', () => {
    const tmp = mkTmp();
    try {
      fs.mkdirSync(path.join(tmp, 'src'));
      fs.writeFileSync(path.join(tmp, 'src', 'foo.ts'), SRC_FOO);
      const stalePath = path.join(tmp, 'src', 'foo-sync.ts');
      fs.writeFileSync(stalePath, '// stale content\n');
      const {code, stderr} = runCli(['--check', '--roots', 'src'], tmp);
      expect(code).toBe(1);
      expect(stderr).toMatch(/out of sync/);
      expect(stderr).toMatch(/foo-sync\.ts/);
      // Original stale content untouched
      expect(fs.readFileSync(stalePath, 'utf8')).toBe('// stale content\n');
    } finally { rmTmp(tmp); }
  });

  test('no codegen sources found — exits 0 with friendly message', () => {
    const tmp = mkTmp();
    try {
      fs.mkdirSync(path.join(tmp, 'src'));
      fs.writeFileSync(path.join(tmp, 'src', 'plain.ts'), `export const x = 1;\n`);
      const {code, stdout} = runCli(['--roots', 'src'], tmp);
      expect(code).toBe(0);
      expect(stdout).toMatch(/no @async-source files found/);
    } finally { rmTmp(tmp); }
  });
});
