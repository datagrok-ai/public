// Async→sync TypeScript codegen via ts-morph AST transform.
// Directive grammar and transformation rules: see README.md.

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
    if (renameMatch) {
      const [, src, dst] = renameMatch;
      if (src === dst)
        throw new Error(`@codegen-rename: '${src}=${dst}' is a no-op (source equals target)`);
      for (const [otherSrc, otherDst] of renames) {
        if (otherDst === dst && otherSrc !== src)
          throw new Error(
            `@codegen-rename: target '${dst}' is already used by ` +
            `'${otherSrc}=${otherDst}' — two sources cannot share a target`);
      }
      renames.set(src, dst);
    }
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
  if (Node.isInterfaceDeclaration(stmt) || Node.isTypeAliasDeclaration(stmt))
    return stmt.getName();
  return null;
}

// Unmatched `@async-only-begin` strips to EOF on purpose: surfaces as a
// downstream parse error rather than silently dropping less than intended.
function stripAsyncOnlyLines(text: string): string {
  const lines = text.split('\n');
  const out: string[] = [];
  let inBlock = false;
  for (const line of lines) {
    if (/\/\/\s*@async-only-begin\b/.test(line)) {
      inBlock = true;
      continue;
    }
    if (/\/\/\s*@async-only-end\b/.test(line)) {
      inBlock = false;
      continue;
    }
    if (inBlock) continue;
    if (/\/\/\s*@async-only\b/.test(line)) continue;
    out.push(line);
  }
  return out.join('\n');
}

function stripDirectiveLines(text: string): string {
  return text
    .split('\n')
    .filter((line) => !/^\s*\/\/\s*@(async-source|codegen-rename)\b/.test(line))
    .join('\n');
}

// One async-modifier removal per call. Re-walks from the SourceFile each
// call: previous edits invalidated descendant references.
function stripOneAsyncModifier(file: SourceFile): boolean {
  let target: Node | undefined;
  file.forEachDescendant((node, traversal) => {
    if (Node.isFunctionDeclaration(node) && node.isAsync()) {
      target = node;
      traversal.stop();
    } else if ((Node.isFunctionExpression(node) || Node.isArrowFunction(node)) && node.isAsync()) {
      target = node;
      traversal.stop();
    }
  });
  if (!target) return false;
  if (Node.isFunctionDeclaration(target)) target.setIsAsync(false);
  else if (Node.isFunctionExpression(target)) target.setIsAsync(false);
  else if (Node.isArrowFunction(target)) target.setIsAsync(false);
  return true;
}

// `replaceWithText` re-parses the file and forgets every descendant node, so
// search from the SourceFile and re-search each pass until nothing matches.
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
    // Routes an async import to its sync sibling: rename the local binding
    // here, the import-rebuild step downstream maps the original name through
    // the rename map. String-literal import names (TS 5+) are skipped.
    for (const imp of file.getImportDeclarations()) {
      for (const ni of imp.getNamedImports()) {
        const aliasNode = ni.getAliasNode();
        const nameNode = ni.getNameNode();
        if (aliasNode) {
          if (aliasNode.getText() === oldName) aliasNode.rename(newName);
        } else if (nameNode.getKind() === SyntaxKind.Identifier) {
          const id = nameNode.asKindOrThrow(SyntaxKind.Identifier);
          if (id.getText() === oldName) id.rename(newName);
        }
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
  namespace?: string;
}

function collectImports(file: SourceFile): ImportInfo[] {
  return file.getImportDeclarations().map((imp) => ({
    moduleSpecifier: imp.getModuleSpecifierValue(),
    named: imp.getNamedImports().map((ni) => ni.getName()),
    namespace: imp.getNamespaceImport()?.getText(),
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
    } else if (Node.isFunctionDeclaration(stmt) || Node.isVariableStatement(stmt) ||
               Node.isInterfaceDeclaration(stmt) || Node.isTypeAliasDeclaration(stmt)) {
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

  // Order: awaits FIRST. After `async` is stripped, `await x` re-parses as
  // the call `await(x)` (since `await` is only a keyword in async functions).
  stripAwaits(work);
  unwrapPromiseTypes(work);
  while (stripOneAsyncModifier(work)) {}

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

  // Async-named imports route through the rename map: the sibling module is
  // expected to export the renamed (sync) name.
  const rebuiltImports: ImportInfo[] = [];
  for (const imp of sourceImports) {
    // Namespace imports aren't subject to `@codegen-rename` — the binding is
    // a local alias.
    if (imp.namespace && referenced.has(imp.namespace) && !declaredInOutput.has(imp.namespace)) {
      rebuiltImports.push({moduleSpecifier: imp.moduleSpecifier, named: [], namespace: imp.namespace});
      continue;
    }
    const used: string[] = [];
    for (const n of imp.named) {
      const resolved = directives.renames.get(n) ?? n;
      if (referenced.has(resolved) && !declaredInOutput.has(resolved))
        used.push(resolved);
    }
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
    `/* eslint-disable */`,
    `// GENERATED — do not edit by hand.`,
    `// Run \`npm run update-codegen\` to regenerate.`,
    `// Source: ./${srcStem}.ts`,
    ``,
  ].join('\n');

  const importLines = rebuiltImports.map((imp) =>
    imp.namespace
      ? `import * as ${imp.namespace} from '${imp.moduleSpecifier}';`
      : `import {${imp.named.join(', ')}} from '${imp.moduleSpecifier}';`,
  ).join('\n');

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
