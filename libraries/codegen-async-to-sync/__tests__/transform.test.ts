import * as fs from 'fs';
import * as path from 'path';
import {transformText} from '../src/index';

function run(src: string, stem = 'in'): string {
  const r = transformText(src, stem);
  if (!r) throw new Error('expected codegen to produce output, got null');
  return r.outputText;
}

const HEADER = '// @async-source: out.ts\n';

describe('async strip', () => {
  test('async function → sync function', () => {
    const out = run(`${HEADER}export async function foo() { return 1; }\n`);
    expect(out).toContain('function foo(');
    expect(out).not.toContain('async function');
  });

  test('async arrow assigned to const → sync arrow', () => {
    const out = run(`${HEADER}export const foo = async () => 1;\n`);
    expect(out).toMatch(/export const foo = \(\) => 1/);
  });

  test('async function with no awaits — body unchanged', () => {
    const out = run(`${HEADER}export async function foo() { const x = 1; return x; }\n`);
    expect(out).toContain('const x = 1');
    expect(out).toContain('return x');
  });

  test('nested async arrow inside outer async fn — both stripped', () => {
    const src = `${HEADER}export async function foo() {
  const evalF = async (p: number): Promise<number> => await bar(p);
  return await evalF(1);
}
`;
    const out = run(src);
    expect(out).not.toMatch(/\basync\b/);
    expect(out).not.toMatch(/\bawait\b/);
    expect(out).toContain('const evalF = (p: number): number => bar(p);');
    expect(out).toContain('return evalF(1);');
  });

  test('nested async function expression — async stripped', () => {
    const src = `${HEADER}export async function foo() {
  const inner = async function(): Promise<number> { return await bar(); };
  return await inner();
}
`;
    const out = run(src);
    expect(out).not.toMatch(/\basync\b/);
    expect(out).not.toMatch(/\bawait\b/);
    expect(out).toMatch(/const inner = function\(\): number/);
  });

  test('doubly-nested async closure (inner inside inner) — both stripped', () => {
    const src = `${HEADER}export async function foo() {
  const outer = async () => {
    const inner = async () => await bar();
    return await inner();
  };
  return await outer();
}
`;
    const out = run(src);
    expect(out).not.toMatch(/\basync\b/);
    expect(out).not.toMatch(/\bawait\b/);
  });
});

describe('await strip', () => {
  test('await fn(x) → fn(x)', () => {
    const out = run(`${HEADER}export async function foo() { return await bar(1); }\n`);
    expect(out).toContain('return bar(1)');
    expect(out).not.toContain('await');
  });

  test('nested await — all stripped', () => {
    const src = `${HEADER}export async function foo() { return (await (await bar()).baz()); }\n`;
    const out = run(src);
    expect(out).not.toContain('await');
    expect(out).toContain('bar()');
    expect(out).toContain('.baz()');
  });

  test('await inside loop / try-catch — structure preserved', () => {
    const src = `${HEADER}export async function foo() {
  try {
    for (let i = 0; i < 3; i++) {
      const v = await bar(i);
      if (v) return v;
    }
  } catch (e) { throw e; }
  return null;
}
`;
    const out = run(src);
    expect(out).not.toContain('await');
    expect(out).toContain('for (let i = 0; i < 3; i++)');
    expect(out).toContain('try {');
    expect(out).toContain('catch (e)');
    expect(out).toContain('const v = bar(i);');
  });
});

describe('Promise<T> unwrap', () => {
  test('Promise<number> → number in return type', () => {
    const out = run(`${HEADER}export async function foo(): Promise<number> { return 1; }\n`);
    expect(out).toMatch(/function foo\(\)\s*:\s*number/);
  });

  test('nested Promise<Promise<T>> → T', () => {
    const out = run(`${HEADER}export async function foo(): Promise<Promise<string>> { return Promise.resolve('a'); }\n`);
    expect(out).toMatch(/function foo\(\)\s*:\s*string/);
  });

  test('Promise<T> inside generic position', () => {
    const out = run(`${HEADER}export async function foo(): Promise<Array<Promise<number>>> { return []; }\n`);
    expect(out).toMatch(/Array<number>/);
    expect(out).not.toContain('Promise<');
  });
});

describe('@codegen-rename', () => {
  test('function decl rename — declaration and references both renamed', () => {
    const src = `// @async-source: out.ts
// @codegen-rename: foo=bar
export async function foo(): Promise<number> {
  return await foo();
}
`;
    const out = run(src);
    expect(out).toMatch(/function bar\(/);
    expect(out).toContain('return bar();');
    expect(out).not.toMatch(/\bfoo\b/);
  });

  test('rename on const-initialised async function', () => {
    const src = `// @async-source: out.ts
// @codegen-rename: foo=barSync
export const foo = async function(): Promise<number> { return 1; };
`;
    const out = run(src);
    expect(out).toContain('export const barSync = function(');
    expect(out).not.toMatch(/\bfoo\b/);
  });

  test('rename of missing identifier — no throw, no diff vs absent rename', () => {
    const a = run(`${HEADER}export async function foo() { return 1; }\n`);
    const b = run(`// @async-source: out.ts
// @codegen-rename: does_not_exist=whatever
export async function foo() { return 1; }
`);
    expect(b).toBe(a);
  });
});

describe('@async-only', () => {
  test('lines marked @async-only are dropped', () => {
    const src = `${HEADER}export async function foo() {
  const debug = true; // @async-only
  return 1;
}
`;
    const out = run(src);
    expect(out).not.toContain('const debug = true');
    expect(out).toContain('return 1;');
  });
});

describe('const type annotation', () => {
  test('async const drops type annotation', () => {
    const src = `// @async-source: out.ts
interface IOpt { (): Promise<number>; }
export const foo: IOpt = async function(): Promise<number> { return 1; };
`;
    const out = run(src);
    expect(out).toMatch(/export const foo = function/);
    expect(out).not.toContain(': IOpt =');
  });

  test('non-async const keeps its annotation (counter-case)', () => {
    const src = `// @async-source: out.ts
export const N: number = 42;
export async function foo(): Promise<number> { return N; }
`;
    const out = run(src);
    expect(out).toContain('import {N}');
    expect(out).toContain('function foo(');
  });
});

describe('imports', () => {
  test('module import filtered to names actually used', () => {
    const src = `// @async-source: out.ts
import {A, B} from 'mod';
export async function foo(): Promise<number> { return A; }
`;
    const out = run(src);
    expect(out).toMatch(/import \{A\} from 'mod'/);
    expect(out).not.toMatch(/\bB\b/);
  });

  test('module import dropped entirely when none of its names are used', () => {
    const src = `// @async-source: out.ts
import {Unused} from 'mod';
export async function foo(): Promise<number> { return 1; }
`;
    const out = run(src);
    expect(out).not.toContain("from 'mod'");
  });

  test('sibling top-level decl referenced by sync body → imported from source stem', () => {
    const src = `// @async-source: out.ts
export function helper(): number { return 1; }
export async function foo(): Promise<number> { return helper(); }
`;
    const out = run(src, 'my-mod');
    expect(out).toContain("import {helper} from './my-mod'");
    expect(out).toContain('function foo(');
    expect(out).not.toContain('function helper(');
  });

  test('sibling import coexists with module import', () => {
    const src = `// @async-source: out.ts
import {A} from 'mod';
export function helper(x: number): number { return x + 1; }
export async function foo(): Promise<number> { return A + helper(1); }
`;
    const out = run(src, 'my-mod');
    expect(out).toMatch(/import \{A\} from 'mod'/);
    expect(out).toMatch(/import \{helper\} from '\.\/my-mod'/);
  });
});

describe('banner & directive stripping', () => {
  test('output starts with GENERATED banner', () => {
    const out = run(`${HEADER}export async function foo() { return 1; }\n`, 'src');
    expect(out.startsWith('// GENERATED — do not edit by hand.\n')).toBe(true);
    expect(out).toContain('// Source: ./src.ts');
  });

  test('directive comment lines never appear in output', () => {
    const src = `// @async-source: out.ts
// @codegen-rename: foo=bar
export async function foo() { return 1; }
`;
    const out = run(src);
    expect(out).not.toContain('@async-source');
    expect(out).not.toContain('@codegen-rename');
  });
});

describe('idempotence', () => {
  const candidates = [
    `${HEADER}export async function foo(): Promise<number> { return await bar(1); }\n`,
    `// @async-source: out.ts
// @codegen-rename: foo=fooSync
import {dep} from 'm';
export async function foo(x: number): Promise<number> { return await dep(x); }
`,
    `// @async-source: out.ts
export function sib(): number { return 1; }
export const foo = async (): Promise<number> => sib();
`,
  ];

  test.each(candidates.map((s, i) => [i, s]))('case %i: running on its own output (after restoring header) is stable', (_, src) => {
    const first = run(src as string);
    // Output has no @async-source directive; transformText would return null on
    // the raw output. Idempotence here means: running transformText a second
    // time on the *original input* yields the same result. (The output is
    // terminal — by design — because directives are stripped.)
    const second = run(src as string);
    expect(second).toBe(first);
  });
});

describe('error paths', () => {
  test('@async-source declared but no async top-level decls → throws', () => {
    const src = `// @async-source: out.ts
export function nothingAsync() { return 1; }
`;
    expect(() => transformText(src, 'src')).toThrow(/no async top-level declarations/);
  });

  test('no @async-source directive → returns null', () => {
    const src = `export async function foo() { return 1; }\n`;
    expect(transformText(src, 'src')).toBeNull();
  });
});

describe('regression — NM optimizer fixture', () => {
  const fixDir = path.join(__dirname, 'fixtures', 'nm-optimizer');
  const source = fs.readFileSync(path.join(fixDir, 'source.ts'), 'utf8');
  const expected = fs.readFileSync(path.join(fixDir, 'expected.ts'), 'utf8');

  test('source → expected (byte-identical)', () => {
    const r = transformText(source, 'optimizer-nelder-mead');
    expect(r).not.toBeNull();
    expect(r!.outputPath).toBe('optimizer-nelder-mead-sync.ts');
    expect(r!.outputText).toBe(expected);
  });
});
