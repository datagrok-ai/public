'use strict';

const USES_U2 = /@datagrok-libraries\/u2|datagrok-api\/u2core|@datagrok-libraries\/ems|DG\.U2\b/;
/** The control class names a rule assumes when it is not given u2's list (check-consumer reads it off the sources). */
const DEFAULT_CONTROLS = /^(Component|Control)$|(Input|Dialog|Form|Menu|Table|List|Tree|Grid|Editor|Picker|Header)$/;
/** u2's documented utility classes (README › Styling): a consumer puts them on its own elements and may style them. */
const PUBLIC_CLASSES = ['u2-adaptive', 'u2-p1', 'u2-p2'];
const PUBLIC_TAIL = `(?!(?:${PUBLIC_CLASSES.map((c) => c.slice(3)).join('|')})(?![\\w-]))`;
/** A u2 class or data-u2 attribute in a selector, the utility classes left out. */
const INTERNAL_SELECTOR = new RegExp(`\\.u2-${PUBLIC_TAIL}|\\[data-u2`);
/** A u2 class name in a class list, the utility classes left out. */
const INTERNAL_CLASS = new RegExp(`(^|\\s)u2-${PUBLIC_TAIL}`);

function libraryOf(source) {
  if (/^@datagrok-libraries\/ems(\/|$)/.test(source))
    return 'ems';
  if (/^@datagrok-libraries\/u2(\/|$)|^datagrok-api\/u2core$/.test(source))
    return 'u2';
  return null;
}

class U2Imports {
  constructor(program) {
    this.names = new Map();
    this.namespaces = new Map();
    for (const st of program.body) {
      if (st.type !== 'ImportDeclaration')
        continue;
      const lib = libraryOf(st.source.value);
      if (lib === null)
        continue;
      for (const s of st.specifiers) {
        const name = s.type === 'ImportSpecifier' ? s.imported.name ?? s.imported.value : 'default';
        if (s.type === 'ImportNamespaceSpecifier')
          this.namespaces.set(s.local.name, lib);
        else
          this.names.set(s.local.name, {lib, name});
      }
    }
  }

  /** `X` imported from u2/ems, `ns.X` over a namespace import of either, or `DG.U2.X`: `{lib, name}`, else null. */
  resolve(node) {
    if (node.type === 'Identifier')
      return this.names.get(node.name) ?? null;
    if (node.type !== 'MemberExpression' || node.computed || node.property.type !== 'Identifier')
      return null;
    const object = node.object;
    if (object.type === 'Identifier' && this.namespaces.has(object.name))
      return {lib: this.namespaces.get(object.name), name: node.property.name};
    if (object.type === 'MemberExpression' && !object.computed && object.property.type === 'Identifier' &&
      object.property.name === 'U2')
      return {lib: 'u2', name: node.property.name};
    return null;
  }

  static usesU2(text) {
    return USES_U2.test(text);
  }
}

const isFunction = (node) => node.type === 'ArrowFunctionExpression' || node.type === 'FunctionExpression' ||
  node.type === 'FunctionDeclaration';

function calleeName(call) {
  const c = call.callee;
  if (c.type === 'Identifier')
    return c.name;
  return c.type === 'MemberExpression' && !c.computed && c.property.type === 'Identifier' ? c.property.name : null;
}

/** `enter` returns false to skip a node's children. */
function walk(root, visitorKeys, enter) {
  const visit = (node) => {
    if (enter(node) === false)
      return;
    for (const key of visitorKeys[node.type] ?? []) {
      const child = node[key];
      for (const c of Array.isArray(child) ? child : [child]) {
        if (c && typeof c.type === 'string')
          visit(c);
      }
    }
  };
  visit(root);
}

/** A literal, a template's fixed parts, or `+` of those. */
function staticText(node) {
  if (node.type === 'Literal' && typeof node.value === 'string')
    return node.value;
  if (node.type === 'TemplateLiteral')
    return node.quasis.map((q) => q.value.cooked ?? q.value.raw).join('');
  if (node.type === 'BinaryExpression' && node.operator === '+')
    return (staticText(node.left) ?? '') + (staticText(node.right) ?? '');
  return null;
}

const sourceOf = (context) => context.sourceCode ?? context.getSourceCode();

/** What the name `id` holds where it is read: the value of the last write to that variable before `id` — its
 * declaration's initializer or an assignment — found through the scopes around `id`; null for a parameter, an
 * import, a name declared nowhere or one not yet written. */
function valueOf(context, id) {
  const source = sourceOf(context);
  for (let scope = source.getScope ? source.getScope(id) : context.getScope(); scope; scope = scope.upper) {
    const variable = scope.set.get(id.name);
    if (variable === undefined)
      continue;
    const writes = variable.references.filter((r) => r.writeExpr && r.writeExpr.range[1] <= id.range[0]);
    return writes.length > 0 ? writes[writes.length - 1].writeExpr : null;
  }
  return null;
}

module.exports = {DEFAULT_CONTROLS, INTERNAL_CLASS, INTERNAL_SELECTOR, PUBLIC_TAIL, U2Imports, calleeName, isFunction,
  sourceOf, staticText, valueOf, walk};
