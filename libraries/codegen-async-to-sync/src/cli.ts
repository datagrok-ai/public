#!/usr/bin/env node
import * as path from 'path';
import {findCodegenSources, processFile, checkFile} from './index';

interface CliArgs {
  roots: string[];
  check: boolean;
}

function parseArgs(argv: string[]): CliArgs {
  const args: CliArgs = {roots: [], check: false};
  let i = 0;
  while (i < argv.length) {
    const a = argv[i];
    if (a === '--check') {
      args.check = true;
      i++;
    } else if (a === '--roots') {
      i++;
      while (i < argv.length && !argv[i].startsWith('--')) {
        args.roots.push(argv[i]);
        i++;
      }
    } else if (a === '--help' || a === '-h') {
      printHelp();
      process.exit(0);
    } else {
      console.error(`codegen-async-to-sync: unknown argument '${a}'`);
      printHelp();
      process.exit(2);
    }
  }
  if (args.roots.length === 0) args.roots = ['src'];
  return args;
}

function printHelp(): void {
  console.log(`Usage: codegen-async-to-sync [--roots <dir>...] [--check]

  --roots <dir>...   Directories to scan for @async-source files (relative to CWD).
                     Default: src
  --check            Verify mode: regenerate in memory and compare with the
                     committed sync sibling. Exit 1 if any file differs. Does not
                     rewrite files.

Files are only processed if their leading comment block contains an
\`@async-source: <output-filename>\` directive.`);
}

function main(): void {
  const args = parseArgs(process.argv.slice(2));
  const cwd = process.cwd();
  const absRoots = args.roots.map((r) => path.resolve(cwd, r));

  const sources = findCodegenSources(absRoots);
  if (sources.length === 0) {
    console.log(`codegen-async-to-sync: no @async-source files found under ${args.roots.join(', ')}`);
    return;
  }

  if (args.check) {
    const drifted: string[] = [];
    for (const src of sources) {
      try {
        const res = checkFile(src);
        if (res?.drift) drifted.push(res.outputPath);
      } catch (err) {
        console.error(`codegen-async-to-sync: ${(err as Error).message}`);
        process.exit(2);
      }
    }
    if (drifted.length > 0) {
      console.error(`codegen-async-to-sync: ${drifted.length} generated file(s) out of sync:`);
      for (const p of drifted) console.error(`  ${path.relative(cwd, p)}`);
      console.error(`Run \`npm run update-codegen\` to regenerate.`);
      process.exit(1);
    }
    console.log(`codegen-async-to-sync: ${sources.length} file(s) up to date.`);
    return;
  }

  for (const src of sources) {
    try {
      const res = processFile(src);
      if (res) console.log(`codegen: wrote ${path.relative(cwd, res.outputPath)}`);
    } catch (err) {
      console.error(`codegen-async-to-sync: ${(err as Error).message}`);
      process.exit(2);
    }
  }
}

main();
