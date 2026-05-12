// Async→sync TypeScript codegen via ts-morph AST transform.
//
// Public surface:
//   transformText(srcText, srcStem) → { outputPath, outputText } | null
//     Pure transform. Returns null if the source has no @async-source directive.
//   processFile(srcPath) → { outputPath, wrote }
//     Reads srcPath, transforms, writes the generated sibling. No-op (returns
//     wrote=false) if the source is not a codegen source.
//   checkFile(srcPath) → { outputPath, drift }
//     Reads srcPath, transforms in memory, compares with what's on disk. Does
//     not write. drift=true if the existing sync file differs from generated.
//   findCodegenSources(roots) → string[]
//     Walks the given root directories and returns absolute paths of files
//     whose leading comment block contains `@async-source`.
//
// Directive grammar (in the source's leading comment block):
//   // @async-source: <output-filename>
//   // @codegen-rename: <old>=<new>          (repeatable)
//   // @async-only                           (end-of-line — drops the host line)

import {
  Project,
  SyntaxKind,
  Node,
  SourceFile,
  Statement,
} from 'ts-morph';
import * as fs from 'fs';
import * as path from 'path';

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

function stripAsyncOnlyLines(text: string): string {
  return text
    .split('\n')
    .filter((line) => !/\/\/\s*@async-only\b/.test(line))
    .join('\n');
}

function stripDirectiveLines(text: string): string {
  return text
    .split('\n')
    .filter((line) => !/^\s*\/\/\s*@(async-source|codegen-rename)\b/.test(line))
    .join('\n');
}

// Find the first async top-level declaration (function or async function-init
// const) and strip its `async` modifier. Returns false when nothing left to
// strip. Re-walks from the SourceFile each call because every previous edit
// has invalidated descendant node references.
function stripOneAsyncModifier(file: SourceFile): boolean {
  for (const stmt of file.getStatements()) {
    if (Node.isFunctionDeclaration(stmt) && stmt.isAsync()) {
      stmt.setIsAsync(false);
      return true;
    }
    if (Node.isVariableStatement(stmt)) {
      for (const decl of stmt.getDeclarations()) {
        const init = decl.getInitializer();
        if (init && (Node.isFunctionExpression(init) || Node.isArrowFunction(init)) && init.isAsync()) {
          init.setIsAsync(false);
          return true;
        }
      }
    }
  }
  return false;
}

// One replacement per pass — `replaceWithText` re-parses the source file, which
// forgets every Node reference except the SourceFile itself. So traverse from
// the SourceFile (stable identity) and re-search each pass until nothing matches.
function stripAwaits(file: SourceFile): void {
  while (true) {
    let target: Node | undefined;
    file.forEachDescendant((node, traversal) => {
      if (node.isKind(SyntaxKind.AwaitExpression)) {
        target = node;
        traversal.stop();
      }
    });
    if (!target) break;
    const expr = target.asKindOrThrow(SyntaxKind.AwaitExpression).getExpression();
    target.replaceWithText(expr.getText());
  }
}

function unwrapPromiseTypes(file: SourceFile): void {
  while (true) {
    let target: Node | undefined;
    file.forEachDescendant((node, traversal) => {
      if (!node.isKind(SyntaxKind.TypeReference)) return;
      const tr = node.asKindOrThrow(SyntaxKind.TypeReference);
      if (tr.getTypeName().getText() !== 'Promise') return;
      if (tr.getTypeArguments().length !== 1) return;
      target = node;
      traversal.stop();
    });
    if (!target) break;
    const tr = target.asKindOrThrow(SyntaxKind.TypeReference);
    target.replaceWithText(tr.getTypeArguments()[0].getText());
  }
}

function applyRenames(file: SourceFile, renames: Map<string, string>): void {
  for (const [oldName, newName] of renames) {
    for (const fn of file.getFunctions()) {
      if (fn.getName() === oldName) fn.rename(newName);
    }
    for (const vs of file.getVariableStatements()) {
      for (const decl of vs.getDeclarations()) {
        if (decl.getName() === oldName) decl.rename(newName);
      }
    }
  }
}

function dropAsyncConstAnnotations(file: SourceFile): void {
  for (const vs of file.getVariableStatements()) {
    for (const decl of vs.getDeclarations()) {
      const init = decl.getInitializer();
      if (!init) continue;
      if (Node.isFunctionExpression(init) || Node.isArrowFunction(init)) {
        if (decl.getTypeNode()) decl.removeType();
      }
    }
  }
}

function collectReferencedNames(scope: Node): Set<string> {
  const out = new Set<string>();
  scope.forEachDescendant((node) => {
    if (node.isKind(SyntaxKind.Identifier)) {
      const parent = node.getParent();
      if (parent && parent.isKind(SyntaxKind.PropertyAccessExpression)) {
        if (parent.asKindOrThrow(SyntaxKind.PropertyAccessExpression).getNameNode() === node)
          return;
      }
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
  named: string[];
}

function collectImports(file: SourceFile): ImportInfo[] {
  return file.getImportDeclarations().map((imp) => ({
    moduleSpecifier: imp.getModuleSpecifierValue(),
    named: imp.getNamedImports().map((ni) => ni.getName()),
  }));
}

export interface TransformResult {
  outputPath: string;
  outputText: string;
}

export function transformText(srcText: string, srcStem: string): TransformResult | null {
  const directives = parseDirectives(srcText);
  if (!directives) return null;

  const cleanedText = stripDirectiveLines(stripAsyncOnlyLines(srcText));

  const project = new Project({
    compilerOptions: {target: 99, module: 99, strict: false, skipLibCheck: true},
    useInMemoryFileSystem: true,
  });
  const work = project.createSourceFile('work.ts', cleanedText);

  const asyncStmts: Statement[] = [];
  const siblingNames = new Set<string>();
  for (const stmt of work.getStatements()) {
    if (isAsyncTopLevel(stmt)) {
      asyncStmts.push(stmt);
    } else if (Node.isFunctionDeclaration(stmt) || Node.isVariableStatement(stmt)) {
      const nm = getDeclaredName(stmt);
      if (nm != null) siblingNames.add(nm);
      if (Node.isVariableStatement(stmt)) {
        for (const d of stmt.getDeclarations()) siblingNames.add(d.getName());
      }
    }
  }
  if (asyncStmts.length === 0)
    throw new Error(`${srcStem}: @async-source declared but no async top-level declarations found`);

  const sourceImports = collectImports(work);

  for (const stmt of [...work.getStatements()]) {
    if (Node.isImportDeclaration(stmt)) continue;
    if (asyncStmts.includes(stmt)) continue;
    stmt.remove();
  }

  // Order matters here. `await` is only a keyword inside async functions; once
  // we strip the `async` modifier, `await x` re-parses as the call expression
  // `await(x)`. So strip awaits FIRST (while the host functions are still
  // async), then unwrap Promise<>, then strip the async modifiers.
  // Each edit re-parses the SourceFile and forgets descendant references, so
  // these helpers traverse from the SourceFile (whose identity survives) and
  // iterate to a fixed point.
  stripAwaits(work);
  unwrapPromiseTypes(work);
  while (stripOneAsyncModifier(work)) { /* iterate to fixed point */ }

  dropAsyncConstAnnotations(work);
  applyRenames(work, directives.renames);

  for (const imp of [...work.getImportDeclarations()]) imp.remove();

  const referenced = new Set<string>();
  for (const stmt of work.getStatements()) {
    for (const name of collectReferencedNames(stmt)) referenced.add(name);
  }

  const declaredInOutput = new Set<string>();
  for (const stmt of work.getStatements()) {
    const nm = getDeclaredName(stmt);
    if (nm) declaredInOutput.add(nm);
  }

  const rebuiltImports: ImportInfo[] = [];
  for (const imp of sourceImports) {
    const used = imp.named.filter((n) => referenced.has(n) && !declaredInOutput.has(n));
    if (used.length === 0) continue;
    rebuiltImports.push({moduleSpecifier: imp.moduleSpecifier, named: used});
  }

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

  const banner = [
    `// GENERATED — do not edit by hand.`,
    `// Run \`npm run update-codegen\` to regenerate.`,
    `// Source: ./${srcStem}.ts`,
    ``,
  ].join('\n');

  const importLines = rebuiltImports.map((imp) =>
    `import {${imp.named.join(', ')}} from '${imp.moduleSpecifier}';`).join('\n');

  const bodyLines = work.getStatements()
    .filter((s) => !Node.isImportDeclaration(s))
    .map((s) => s.getText())
    .join('\n\n');

  const outputText = `${banner}${importLines}\n\n${bodyLines}\n`;
  return {outputPath: directives.outputPath, outputText};
}

export interface ProcessFileResult {
  outputPath: string;
  wrote: boolean;
}

export function processFile(srcPath: string): ProcessFileResult | null {
  const srcText = fs.readFileSync(srcPath, 'utf8');
  const srcStem = path.basename(srcPath, '.ts');
  const result = transformText(srcText, srcStem);
  if (!result) return null;
  const outAbs = path.join(path.dirname(srcPath), result.outputPath);
  fs.writeFileSync(outAbs, result.outputText, 'utf8');
  return {outputPath: outAbs, wrote: true};
}

export interface CheckFileResult {
  outputPath: string;
  drift: boolean;
  expected: string;
  actual: string;
}

export function checkFile(srcPath: string): CheckFileResult | null {
  const srcText = fs.readFileSync(srcPath, 'utf8');
  const srcStem = path.basename(srcPath, '.ts');
  const result = transformText(srcText, srcStem);
  if (!result) return null;
  const outAbs = path.join(path.dirname(srcPath), result.outputPath);
  const existing = fs.existsSync(outAbs) ? fs.readFileSync(outAbs, 'utf8') : '';
  return {
    outputPath: outAbs,
    drift: existing !== result.outputText,
    expected: result.outputText,
    actual: existing,
  };
}

export function findCodegenSources(roots: string[]): string[] {
  const out: string[] = [];
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
    const head = fs.readFileSync(full, 'utf8').slice(0, 512);
    if (head.includes('@async-source')) acc.push(full);
  }
}
