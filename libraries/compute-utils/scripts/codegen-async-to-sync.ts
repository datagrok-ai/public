// Psycopg-style async→sync codegen.
//
// Walks files in this library that carry a `// @async-source:` directive,
// derives a sync mirror from each, and writes it next to the source. Run
// manually via `npm run update-codegen`; CI verifies no drift via
// `npm run codegen:check`.
//
// Recognised directives (in the source's leading comment block):
//   // @async-source: <output-filename>    declares this file as a codegen
//                                          input and names the emitted sibling
//   // @codegen-rename: <old>=<new>        renames a top-level identifier in
//                                          the sync output (repeatable)
//   // @async-only                         end-of-line marker — drops the
//                                          host line from the sync output
//
// Transformation rules (deterministic, mechanical):
//   - `async function` / `async` arrow → strip the `async`
//   - `await expr`                     → `expr`
//   - `Promise<T>`                     → `T`
//   - Apply @codegen-rename to declared name and all references
//   - Top-level declarations that aren't async → not copied to sync output;
//     if the copied code references them they're imported from the source.
//
// Single-source-of-truth invariants enforced by `codegen:check` in CI.

import {
  Project,
  SyntaxKind,
  Node,
  SourceFile,
  Statement,
  FunctionDeclaration,
  FunctionExpression,
  ArrowFunction,
} from 'ts-morph';
import * as fs from 'fs';
import * as path from 'path';

const ROOT = path.resolve(__dirname, '..');

interface Directives {
  outputPath: string;
  renames: Map<string, string>;
}

function parseDirectives(sourceText: string): Directives | null {
  const lines = sourceText.split('\n');
  const renames = new Map<string, string>();
  let outputPath: string | null = null;
  for (const raw of lines) {
    const line = raw.trim();
    if (line === '') continue;
    if (!line.startsWith('//')) break;
    const asyncSource = /@async-source:\s*(\S+)/.exec(line);
    if (asyncSource) outputPath = asyncSource[1];
    const renameMatch = /@codegen-rename:\s*([A-Za-z_$][\w$]*)\s*=\s*([A-Za-z_$][\w$]*)/.exec(line);
    if (renameMatch) renames.set(renameMatch[1], renameMatch[2]);
  }
  return outputPath ? {outputPath, renames} : null;
}

// "Async" top-level statements get copied (with async/await stripped) into
// the sync output. Everything else is left in the source file and imported
// back from there when the sync code references it.
function isAsyncTopLevel(stmt: Statement): boolean {
  if (Node.isFunctionDeclaration(stmt))
    return stmt.isAsync();
  if (Node.isVariableStatement(stmt)) {
    for (const decl of stmt.getDeclarations()) {
      const init = decl.getInitializer();
      if (!init) continue;
      if ((Node.isFunctionExpression(init) || Node.isArrowFunction(init)) && init.isAsync())
        return true;
    }
  }
  return false;
}

function getDeclaredName(stmt: Statement): string | null {
  if (Node.isFunctionDeclaration(stmt)) return stmt.getName() ?? null;
  if (Node.isVariableStatement(stmt)) {
    const decls = stmt.getDeclarations();
    return decls.length === 1 ? decls[0].getName() : null;
  }
  return null;
}

// Strip @async-only lines (after their position survives transforms is hard
// to track via AST). Source-level pre-processing of text is the simplest
// honest approach for this directive — psycopg does the same trick.
function stripAsyncOnlyLines(text: string): string {
  return text
    .split('\n')
    .filter((line) => !/\/\/\s*@async-only\b/.test(line))
    .join('\n');
}

// Strip the directive comment lines from the cloned source so they don't
// re-trigger codegen if someone copies the output.
function stripDirectiveLines(text: string): string {
  return text
    .split('\n')
    .filter((line) => !/^\s*\/\/\s*@(async-source|codegen-rename)\b/.test(line))
    .join('\n');
}

// Strip `async` from a function decl / function expression / arrow.
function stripAsync(fn: FunctionDeclaration | FunctionExpression | ArrowFunction): void {
  if (!fn.isAsync()) return;
  fn.setIsAsync(false);
}

// Replace every `await expr` with just `expr`.
function stripAwaits(scope: Node): void {
  // Walk descendants in reverse so replacing children doesn't invalidate
  // ancestor iteration. forEachDescendant doesn't expose this; collect first.
  const awaits: Node[] = [];
  scope.forEachDescendant((node) => {
    if (node.isKind(SyntaxKind.AwaitExpression)) awaits.push(node);
  });
  for (const aw of awaits) {
    const expr = aw.asKindOrThrow(SyntaxKind.AwaitExpression).getExpression();
    aw.replaceWithText(expr.getText());
  }
}

// Replace every `Promise<T>` type reference with `T`. Handles nested
// Promise<Promise<...>> by iterating to a fixed point.
function unwrapPromiseTypes(scope: Node): void {
  while (true) {
    const refs: Node[] = [];
    scope.forEachDescendant((node) => {
      if (node.isKind(SyntaxKind.TypeReference) &&
        node.asKindOrThrow(SyntaxKind.TypeReference).getTypeName().getText() === 'Promise') {
        refs.push(node);
      }
    });
    if (refs.length === 0) break;
    for (const ref of refs) {
      const tr = ref.asKindOrThrow(SyntaxKind.TypeReference);
      const args = tr.getTypeArguments();
      if (args.length !== 1) continue;
      ref.replaceWithText(args[0].getText());
    }
  }
}

// Rename top-level identifiers and any reference to them in the cloned code.
function applyRenames(file: SourceFile, renames: Map<string, string>): void {
  for (const [oldName, newName] of renames) {
    // Function declarations
    for (const fn of file.getFunctions()) {
      if (fn.getName() === oldName) fn.rename(newName);
    }
    // Top-level variable declarations
    for (const vs of file.getVariableStatements()) {
      for (const decl of vs.getDeclarations()) {
        if (decl.getName() === oldName) decl.rename(newName);
      }
    }
  }
}

// Drop the const-level type annotation when the value is an async function.
// Example: `export const optimizeNM: IOptimizer = async function(...)` →
//          `export const optimizeNMSync = function(...)`
// The function's own return-type annotation survives intact.
function dropAsyncConstAnnotations(file: SourceFile): void {
  for (const vs of file.getVariableStatements()) {
    for (const decl of vs.getDeclarations()) {
      const init = decl.getInitializer();
      if (!init) continue;
      // After stripAsync(), the init is no longer flagged async — so we
      // check the *original* shape by inspecting init kind.
      if (Node.isFunctionExpression(init) || Node.isArrowFunction(init)) {
        if (decl.getTypeNode()) decl.removeType();
      }
    }
  }
}

// What identifiers does this scope reference (free variables)?
// We don't need rigorous scope analysis — we just collect any Identifier
// at use-position and let the caller cross-check against known top-level
// names. Accidental over-collection (e.g. property names) is harmless;
// they won't match top-level declarations.
function collectReferencedNames(scope: Node): Set<string> {
  const out = new Set<string>();
  scope.forEachDescendant((node) => {
    if (node.isKind(SyntaxKind.Identifier)) {
      const parent = node.getParent();
      // Skip property-access right-hand sides (e.g. `obj.foo` → foo isn't a free var)
      if (parent && parent.isKind(SyntaxKind.PropertyAccessExpression)) {
        if (parent.asKindOrThrow(SyntaxKind.PropertyAccessExpression).getNameNode() === node)
          return;
      }
      // Skip property assignment / signature property names
      if (parent && (
        parent.isKind(SyntaxKind.PropertyAssignment) ||
        parent.isKind(SyntaxKind.PropertySignature) ||
        parent.isKind(SyntaxKind.ShorthandPropertyAssignment)
      )) {
        const nameNode = (parent as any).getNameNode?.();
        if (nameNode === node) return;
      }
      out.add(node.getText());
    }
  });
  return out;
}

interface ImportInfo {
  moduleSpecifier: string;
  named: string[]; // names imported (after possible renames-as-import)
}

function collectImports(file: SourceFile): ImportInfo[] {
  return file.getImportDeclarations().map((imp) => ({
    moduleSpecifier: imp.getModuleSpecifierValue(),
    named: imp.getNamedImports().map((ni) => ni.getName()),
  }));
}

function processFile(srcPath: string): void {
  const srcText = fs.readFileSync(srcPath, 'utf8');
  const directives = parseDirectives(srcText);
  if (!directives) return; // not a codegen source

  const outRel = directives.outputPath;
  const outPath = path.join(path.dirname(srcPath), outRel);
  const srcStem = path.basename(srcPath, '.ts');

  // Strip @async-only lines BEFORE parsing — keeps the AST stable so we
  // don't have to re-find positions after edits.
  const cleanedText = stripDirectiveLines(stripAsyncOnlyLines(srcText));

  // Fresh in-memory project so the codegen never modifies the on-disk source.
  const project = new Project({
    compilerOptions: {target: 99, module: 99, strict: false, skipLibCheck: true},
    useInMemoryFileSystem: true,
  });
  const work = project.createSourceFile('work.ts', cleanedText);

  // Identify async top-level statements; the rest is "siblings" referenced
  // back from the source file in the sync output's imports.
  const asyncStmts: Statement[] = [];
  const siblingNames = new Set<string>(); // names declared in source that we skip
  for (const stmt of work.getStatements()) {
    if (isAsyncTopLevel(stmt)) {
      asyncStmts.push(stmt);
    } else if (Node.isFunctionDeclaration(stmt) || Node.isVariableStatement(stmt)) {
      const nm = getDeclaredName(stmt);
      if (nm != null) siblingNames.add(nm);
      // Multi-declarator VariableStatement: capture each name.
      if (Node.isVariableStatement(stmt)) {
        for (const d of stmt.getDeclarations()) siblingNames.add(d.getName());
      }
    }
  }
  if (asyncStmts.length === 0) {
    throw new Error(`${srcPath}: declared @async-source but no async top-level declarations found`);
  }

  // Original imports (before we drop statements)
  const sourceImports = collectImports(work);

  // Remove non-async top-level statements so only the async ones survive
  // in the working source file.
  for (const stmt of [...work.getStatements()]) {
    if (Node.isImportDeclaration(stmt)) continue;
    if (asyncStmts.includes(stmt)) continue;
    stmt.remove();
  }

  // Drop async modifiers and await expressions.
  for (const stmt of work.getStatements()) {
    if (Node.isFunctionDeclaration(stmt)) {
      stripAsync(stmt);
      stripAwaits(stmt);
      unwrapPromiseTypes(stmt);
    } else if (Node.isVariableStatement(stmt)) {
      for (const decl of stmt.getDeclarations()) {
        const init = decl.getInitializer();
        if (init && (Node.isFunctionExpression(init) || Node.isArrowFunction(init))) {
          stripAsync(init);
          stripAwaits(init);
          unwrapPromiseTypes(init);
        }
      }
    }
  }

  // Drop the const-level type annotation when the value was async (e.g. IOptimizer).
  dropAsyncConstAnnotations(work);

  // Apply renames AFTER async-strip so the renamer sees a stable AST.
  applyRenames(work, directives.renames);

  // Drop the original imports — we'll rebuild them based on what the sync
  // body actually references.
  for (const imp of [...work.getImportDeclarations()]) imp.remove();

  // What does the sync body actually reference now?
  const referenced = new Set<string>();
  for (const stmt of work.getStatements()) {
    for (const name of collectReferencedNames(stmt)) referenced.add(name);
  }

  // Identifiers declared inside the sync output itself (so we don't try to
  // re-import them from elsewhere).
  const declaredInOutput = new Set<string>();
  for (const stmt of work.getStatements()) {
    const nm = getDeclaredName(stmt);
    if (nm) declaredInOutput.add(nm);
  }

  // For renamed identifiers, the *new* name is what appears in `referenced`.
  // Map them back to detect originating modules. (We don't need this for
  // siblings because renaming only applies to top-level decls we keep.)

  // 1. Rebuild original-module imports, filtered to names the sync body uses.
  const rebuiltImports: ImportInfo[] = [];
  for (const imp of sourceImports) {
    const used = imp.named.filter((n) => referenced.has(n) && !declaredInOutput.has(n));
    if (used.length === 0) continue;
    rebuiltImports.push({moduleSpecifier: imp.moduleSpecifier, named: used});
  }

  // 2. Sibling import: anything referenced that's a top-level name in the
  //    source file that we *skipped*.
  const siblingsImported: string[] = [];
  for (const name of referenced) {
    if (declaredInOutput.has(name)) continue;
    if (siblingNames.has(name)) siblingsImported.push(name);
  }
  if (siblingsImported.length > 0) {
    rebuiltImports.push({
      moduleSpecifier: `./${srcStem}`,
      named: siblingsImported.sort(),
    });
  }

  // Render the final output.
  const banner = [
    `// GENERATED — do not edit by hand.`,
    `// Run \`npm run update-codegen\` (in libraries/compute-utils) to regenerate.`,
    `// Source: ./${srcStem}.ts`,
    ``,
  ].join('\n');

  const importLines = rebuiltImports.map((imp) =>
    `import {${imp.named.join(', ')}} from '${imp.moduleSpecifier}';`).join('\n');

  const bodyLines = work.getStatements()
    .filter((s) => !Node.isImportDeclaration(s))
    .map((s) => s.getText())
    .join('\n\n');

  const finalText = `${banner}${importLines}\n\n${bodyLines}\n`;
  fs.writeFileSync(outPath, finalText, 'utf8');
  console.log(`codegen: wrote ${path.relative(ROOT, outPath)}`);
}

function findCodegenSources(): string[] {
  const out: string[] = [];
  // Limit the walk to source dirs we care about; skip node_modules etc.
  const roots = [path.join(ROOT, 'function-views'), path.join(ROOT, 'webworkers'),
    path.join(ROOT, 'shared-utils'), path.join(ROOT, 'shared-components')];
  for (const root of roots) {
    if (!fs.existsSync(root)) continue;
    walk(root, out);
  }
  return out;
}

function walk(dir: string, acc: string[]): void {
  for (const entry of fs.readdirSync(dir, {withFileTypes: true})) {
    const full = path.join(dir, entry.name);
    if (entry.isDirectory()) {
      if (entry.name === 'node_modules' || entry.name === 'dist') continue;
      walk(full, acc);
      continue;
    }
    if (!entry.isFile()) continue;
    if (!entry.name.endsWith('.ts')) continue;
    // Cheap header sniff — full directive parse only after this passes.
    const head = fs.readFileSync(full, 'utf8').slice(0, 512);
    if (head.includes('@async-source')) acc.push(full);
  }
}

function main(): void {
  const sources = findCodegenSources();
  if (sources.length === 0) {
    console.log('codegen: no @async-source files found.');
    return;
  }
  for (const s of sources) processFile(s);
}

main();
