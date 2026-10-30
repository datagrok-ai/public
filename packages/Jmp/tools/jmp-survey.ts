// Surveys the saved analyses in a folder of .jmp files: which JMP platforms they launch, which become
// Datagrok viewers, what doesn't come across, and the JSL vocabulary each platform uses (JMP_PLUGIN.md §11).
// Made to run where the files are, by people without the package:
//
//   node jmp-survey.mjs <folder> [--summary-only]
//
// --summary-only prints counts and JSL keywords only: no file paths, script names, column names or values,
// so the output can be shared outside the organization that owns the files.
//
// Bundle: npx esbuild tools/jmp-survey.ts --bundle --platform=node --format=esm --outfile=jmp-survey.mjs

import {readdirSync, statSync} from 'node:fs';
import {readFile} from 'node:fs/promises';
import {join, relative} from 'node:path';

import {readJmp} from '../src/formats/jmp/jmp-parser';
import {JslNode, jslKey} from '../src/formats/jmp/jsl';
import {isSourceScript} from '../src/formats/jmp/jmp-scripts';
import {convertScripts} from '../src/formats/jmp/jsl-to-viewers';

const MAX_DEPTH = 4;
const SEP = ' › ';

class Counter {
  readonly counts = new Map<string, number>();
  add(key: string, n = 1): void {
    this.counts.set(key, (this.counts.get(key) ?? 0) + n);
  }
  sorted(): [string, number][] {
    return [...this.counts].sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]));
  }
}

function walk(dir: string): string[] {
  return readdirSync(dir).flatMap((name) => {
    const path = join(dir, name);
    return statSync(path).isDirectory() ? walk(path) : name.toLowerCase().endsWith('.jmp') ? [path] : [];
  });
}

/** The platform's call: the outermost call with the platform's name. */
function platformCall(node: JslNode, key: string): Extract<JslNode, {type: 'call'}> | undefined {
  if (node.type === 'call' && jslKey(node.name) === key)
    return node;
  for (const child of children(node)) {
    const found = platformCall(child, key);
    if (found)
      return found;
  }
  return undefined;
}

function children(node: JslNode): JslNode[] {
  switch (node.type) {
  case 'call': return node.args;
  case 'list': return node.items;
  case 'index': return [node.target, ...node.args];
  case 'unary': return [node.operand];
  case 'binary': return [node.left, node.right];
  default: return [];
  }
}

/** Every keyword path under the platform's call, e.g. "Graph Builder › Elements › Points". Calls and bare
 * names are JSL keywords; strings, numbers and column references are the user's and are left out. */
function keywordPaths(root: Extract<JslNode, {type: 'call'}>, spelling: Map<string, string>): Set<string> {
  const out = new Set<string>();
  const visit = (node: JslNode, path: string[]): void => {
    if ((node.type === 'call' || node.type === 'name') && path.length < MAX_DEPTH) {
      const key = jslKey(node.name);
      if (!spelling.has(key))
        spelling.set(key, node.name);
      const next = [...path, key];
      out.add(next.join(SEP));
      if (node.type === 'call')
        node.args.forEach((a) => visit(a, next));
      return;
    }
    children(node).forEach((c) => visit(c, path));
  };
  visit(root, []);
  return out;
}

/** Takes the user's words out of a message: anything quoted longer than a character. The converter quotes
 * every column name it writes, so JSL's own words, such as "Fit Model", are left alone even when a column
 * shares one. */
function redact(message: string): string {
  return message.replace(/"[^"]{2,}"/g, '"…"');
}

async function main(): Promise<void> {
  const args = process.argv.slice(2);
  const summaryOnly = args.includes('--summary-only');
  const dir = args.find((a) => !a.startsWith('--'));
  if (!dir) {
    console.error('usage: node jmp-survey.mjs <folder> [--summary-only]');
    process.exit(2);
  }

  const files = walk(dir);
  let read = 0; let withAnalyses = 0;
  const failures = new Counter();
  const versions = new Counter();
  const platforms = new Map<string, {name: string, scripts: number, converted: number, unread: number}>();
  const reasons = new Counter();
  const notes = new Map<string, Counter>();
  const vocabulary = new Map<string, Counter>();
  const spelling = new Map<string, string>();
  let scripts = 0; let unreadable = 0; let noPlatform = 0;

  for (const file of files) {
    const shown = relative(dir, file);
    try {
      const table = await readJmp(new Uint8Array(await readFile(file)));
      read++;
      versions.add(table.info.jmpVersion.split('.')[0] || '?');
      const saved = table.scripts.filter((s) => !isSourceScript(s));
      const analyses = convertScripts(saved, table.columns.filter((c) => !c.decodeError && c.kind !== 'rowstate')
        .map((c) => ({name: c.info.name, numeric: c.kind !== 'string'})));
      if (analyses.length > 0)
        withAnalyses++;
      analyses.forEach((a, i) => {
        scripts++;
        const script = saved[i];
        if (script.error) {
          unreadable++;
          reasons.add(redact(a.reason!));
        }
        if (!a.platform) {
          if (!script.error) {
            noPlatform++;
            reasons.add(redact(a.reason!));
          }
          return;
        }
        const key = jslKey(a.platform);
        const p = platforms.get(key) ?? {name: a.platform, scripts: 0, converted: 0, unread: 0};
        p.scripts++;
        if (a.viewers.length > 0)
          p.converted++;
        else
          reasons.add(`${a.platform}: ${redact(a.reason!)}`);
        platforms.set(key, p);
        const n = notes.get(key) ?? new Counter();
        for (const note of new Set(a.notes.map((x) => redact(x))))
          n.add(note);
        notes.set(key, n);
        const root = script.body && platformCall(script.body, key);
        if (root) {
          const v = vocabulary.get(key) ?? new Counter();
          for (const path of keywordPaths(root, spelling))
            v.add(path);
          vocabulary.set(key, v);
        }
        if (!summaryOnly) {
          console.log(`${shown}: ${a.name} [${a.platform}] ` +
            (a.viewers.length > 0 ? a.viewers.map((x) => x.type).join(' + ') : `not converted: ${a.reason}`));
        }
      });
    } catch (e) {
      const message = e instanceof Error ? e.message : String(e);
      failures.add(redact(message));
      if (!summaryOnly)
        console.log(`${shown}: FAILED ${message}`);
    }
  }

  const show = (path: string): string => path.split(SEP).map((k) => spelling.get(k) ?? k).join(SEP);
  const indent = (path: string): string => '  '.repeat(path.split(SEP).length - 1);

  console.log(`\n== Files\n${files.length} .jmp files: ${read} read, ${files.length - read} failed, ` +
    `${withAnalyses} with saved analyses`);
  for (const [m, n] of failures.sorted())
    console.log(`  ${n} × ${m}`);
  console.log('Saved by JMP version: ' + versions.sorted().map(([v, n]) => `${v}: ${n}`).join(', '));

  console.log(`\n== Saved analyses\n${scripts} scripts (the "Source" re-import left out): ` +
    `${unreadable} unreadable, ${noPlatform} launch no JMP platform`);
  const byUse = [...platforms.values()].sort((a, b) => b.scripts - a.scripts);
  for (const p of byUse)
    console.log(`  ${p.name}: ${p.scripts}, ${p.converted} become Datagrok viewers`);
  if (reasons.counts.size > 0) {
    console.log('Not converted, by reason:');
    for (const [m, n] of reasons.sorted())
      console.log(`  ${n} × ${m}`);
  }

  console.log('\n== What doesn\'t come across, by platform (scripts with each note)');
  for (const p of byUse) {
    const n = notes.get(jslKey(p.name));
    if (!n || n.counts.size === 0)
      continue;
    console.log(p.name);
    for (const [m, c] of n.sorted())
      console.log(`  ${c} × ${m}`);
  }

  console.log(`\n== JSL vocabulary, by platform (scripts using each keyword, to depth ${MAX_DEPTH})`);
  for (const p of byUse) {
    const v = vocabulary.get(jslKey(p.name));
    if (!v)
      continue;
    // as a tree: each path under its parent, siblings by use
    const paths = v.sorted();
    const depth = (path: string): number => path.split(SEP).length;
    const childrenOf = (parent: string): [string, number][] =>
      paths.filter(([path]) => path.startsWith(parent + SEP) && depth(path) === depth(parent) + 1);
    const print = (path: string, n: number): void => {
      console.log(`${indent(path)}${show(path).split(SEP).pop()}: ${n}`);
      childrenOf(path).forEach(([c, m]) => print(c, m));
    };
    paths.filter(([path]) => !path.includes(SEP)).forEach(([path, n]) => print(path, n));
  }
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
