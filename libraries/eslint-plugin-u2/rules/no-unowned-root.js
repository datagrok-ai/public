'use strict';

const {DEFAULT_CONTROLS, U2Imports, calleeName, isFunction, sourceOf, valueOf, walk} = require('../lib/u2-imports');

const isScope = (node) => node.type === 'Identifier' && node.name === 'Scope' ||
  node.type === 'MemberExpression' && !node.computed && node.property.type === 'Identifier' &&
  node.property.name === 'Scope';

const isNothing = (node) => node === undefined ||
  node.type === 'Identifier' && node.name === 'undefined' ||
  node.type === 'Literal' && node.value === null && node.raw === 'null' ||
  node.type === 'UnaryExpression' && node.operator === 'void';

// what an owner verb is handed runs with the owner ambient (`dialog.add(() => …)` is a dialog's builder)
const OWNER_VERBS = new Set(['build', 'runInScope', 'runWith', 'task', 'effect', 'add']);
// these run their callback before they answer, inside whoever called them
const SYNC_CALLBACKS = new Set(['map', 'flatMap', 'forEach', 'filter', 'reduce', 'find', 'some', 'every', 'sort']);
const WRAPPERS = new Set(['TSAsExpression', 'TSNonNullExpression', 'TSSatisfiesExpression', 'TSTypeAssertion',
  'ChainExpression']);
// how far a helper's callers are followed back
const DEPTH = 4;
const VIEW_TYPE = /\b(?:View|ViewBase|TableView)\b/;
// what makes or shows a platform view: `DG.View.fromRoot(…)`, `DG.View.create()`, `grok.shell.addView(…)`
const VIEW_CALLS = new Set(['fromRoot', 'addView', 'addPreview', 'addTablePreview', 'addTableView']);

const keyName = (key) => key.type === 'Identifier' ? key.name : key.type === 'Literal' ? String(key.value) : null;

const unwrap = (node) => {
  let e = node;
  while (e && (WRAPPERS.has(e.type) || e.type === 'AwaitExpression'))
    e = e.type === 'AwaitExpression' ? e.argument : e.expression;
  return e;
};

const nameOf = (node) => node.type === 'Identifier' ? node.name :
  node.type === 'MemberExpression' && !node.computed ? keyName(node.property) : null;

/** What `expr` hands on as a value: itself, or the parts of a literal, a ternary or an `a ?? b` it is made of. */
function valueParts(expr, out = []) {
  if (!expr)
    return out;
  if (WRAPPERS.has(expr.type))
    return valueParts(expr.expression, out);
  if (expr.type === 'ConditionalExpression') {
    valueParts(expr.consequent, out);
    valueParts(expr.alternate, out);
  }
  else if (expr.type === 'LogicalExpression') {
    valueParts(expr.left, out);
    valueParts(expr.right, out);
  }
  else if (expr.type === 'ArrayExpression') {
    for (const e of expr.elements)
      valueParts(e?.type === 'SpreadElement' ? e.argument : e, out);
  }
  else if (expr.type === 'ObjectExpression') {
    for (const p of expr.properties)
      valueParts(p.type === 'Property' ? p.value : p.argument, out);
  }
  else
    out.push(expr);
  return out;
}

module.exports = {
  meta: {
    type: 'problem',
    docs: {description: 'A new ownership tree is started only on purpose, as a named root, and a control is built ' +
      'under its owner'},
    messages: {
      runWith: '`Scope.runWith(undefined, …)` builds a tree nobody owns, which never ends. Build it under its ' +
        'owner (the view or control that shows it), or start a named root on purpose: `Scope.root(\'<what>\')`.',
      newScope: '`new Scope()` without a parent starts a tree nobody owns, which never ends. Give it its owner\'s ' +
        'scope as the parent, or start a named root on purpose: `Scope.root(\'<what>\')`.',
      unownedControl: '`{{what}}` is built with no owner: what it holds (an observer, a subscription) is never ' +
        'released, and u2 strict mode throws on it. Build it under the view or control that shows it — ' +
        '`ViewHost.of(view).build(() => {{what}})`, `owner.build(() => …)`, inside `component((self) => …)` — or ' +
        'under a named root on purpose: `Scope.root(\'<what>\').build(…)`.',
      unownedCall: '`{{what}}` is built with no owner: what it holds (an observer, a subscription, a timer) is ' +
        'never released. Build it under the view or control that shows it — `ViewHost.create({name, content: () => ' +
        '{{what}}})`, `ViewHost.of(view).build(() => {{what}})`, inside `component((self) => …)` — or under a named ' +
        'root on purpose: `Scope.root(\'<what>\').build(…)`.',
      hostView: '`{{what}}` made into a platform view is never disposed: the view closes and what the control holds ' +
        'lives on. Make the view with `ViewHost.create({name, content: () => …})` and show `host.view` ' +
        '(`grok.shell.addView(host.view)`, or answer it from the app function).',
      hostWidget: '`{{what}}` handed to `DG.Widget.fromRoot` is never disposed when the panel closes. Hand over ' +
        '`host(control)` from u2/dg, which the platform\'s close disposes: `DG.Widget.fromRoot(host(control))`.',
      hostDialog: '`{{what}}` added to a platform dialog is never disposed when the dialog closes. Use u2\'s ' +
        '`Dialog` (`new Dialog(title).add(() => …).show()`), or add `host(control)` from u2/dg.',
      hostDock: '`{{what}}` docked by hand outlives its panel and its view. Dock with ' +
        '`ViewHost.of(view).dock(control, place, {title})`: closing the panel or the view disposes the control.',
      afterAwait: '`{{what}}` is built after an `await`, where nothing is building, so it has no owner. Build it ' +
        'back in its owner: `owner.task(async (t) => { await …; t.build(() => {{what}}); })`.',
    },
    schema: [{type: 'object', properties: {
      controls: {type: 'array', items: {type: 'string'}},
      roots: {type: 'array', items: {type: 'string'}},
      factories: {type: 'array', items: {type: 'string'}},
    }, additionalProperties: false}],
  },
  create(context) {
    const options = context.options[0] ?? {};
    const controls = options.controls ? new Set(options.controls) : null;
    const roots = new Set(options.roots ?? ['Dialog', 'ViewHost']);
    const factories = new Set(options.factories ?? []);
    const source = sourceOf(context);
    const keys = source.visitorKeys;
    let imports;
    const built = [];
    const calls = [];
    const handed = [];
    // this file's functions that answer a control they build: `const panel = () => component(…)`
    const localBuilders = new Set();

    const isControl = (name) => !roots.has(name) &&
      (controls === null ? DEFAULT_CONTROLS.test(name) : controls.has(name));

    /** `new X(…)` of a u2 or ems control, a call of a control factory imported from either or of a local builder;
     * else null. */
    const builtControl = (node) => {
      const ref = imports.resolve(node.callee);
      if (ref === null) {
        return node.type === 'CallExpression' && node.callee.type === 'Identifier' &&
          localBuilders.has(node.callee.name) ? `${node.callee.name}(…)` : null;
      }
      if (node.type === 'NewExpression')
        return isControl(ref.name) ? `new ${ref.name}(…)` : null;
      return factories.has(ref.name) ? `${ref.name}(…)` : null;
    };

    // runs what it is handed with an owner ambient: an owner verb, or a u2 or ems API (`component(…)`,
    // `ViewHost.create({content})`, `new List({render})`), which builds it under its own scope
    const owns = (call) => {
      if (call.type === 'CallExpression' && OWNER_VERBS.has(calleeName(call)))
        return true;
      const c = call.callee;
      return imports.resolve(c) !== null || c.type === 'MemberExpression' && imports.resolve(c.object) !== null;
    };

    const runsNow = (call) => call.type === 'CallExpression' && SYNC_CALLBACKS.has(calleeName(call));

    const controlClass = (cls) => cls.superClass != null && (imports.resolve(cls.superClass) !== null ||
      isControl(nameOf(cls.superClass) ?? ''));

    /** The call or `new` that `fn` is handed to, directly or inside an options object or an array. */
    const handedTo = (fn) => {
      let arg = fn;
      while (arg.parent.type === 'Property' && arg.parent.value === arg || arg.parent.type === 'ObjectExpression' ||
        arg.parent.type === 'ArrayExpression')
        arg = arg.parent;
      const p = arg.parent;
      return (p.type === 'CallExpression' || p.type === 'NewExpression') && p.arguments.includes(arg) ? p : null;
    };

    const afterAwait = (fn, node) => {
      if (!fn.async)
        return false;
      let found = false;
      walk(fn.body, keys, (n) => {
        if (found || n !== fn.body && isFunction(n))
          return false;
        if (n.type === 'AwaitExpression' && n.range[1] <= node.range[0])
          found = true;
        return !found;
      });
      return found;
    };

    /** What `fn` answers: its expression body, or the arguments of its own `return`s. */
    const answered = (fn) => {
      if (fn.body.type !== 'BlockStatement')
        return [fn.body];
      const out = [];
      walk(fn.body, keys, (n) => {
        if (isFunction(n))
          return false;
        if (n.type === 'ReturnStatement' && n.argument)
          out.push(n.argument);
      });
      return out;
    };

    /** Whether `fn` answers what `node` builds — it is a factory, and its caller owns the result. */
    const returns = (fn, node) => {
      const values = answered(fn).flatMap((v) => valueParts(v));
      if (values.includes(node))
        return true;
      let at = node;
      while (WRAPPERS.has(at.parent.type))
        at = at.parent;
      const decl = at.parent;
      return decl.type === 'VariableDeclarator' && decl.init === at && decl.id.type === 'Identifier' &&
        values.some((v) => v.type === 'Identifier' && v.name === decl.id.name);
    };

    const makesView = (expr) => {
      let e = expr;
      while (e && WRAPPERS.has(e.type))
        e = e.expression;
      if (e?.type === 'AwaitExpression')
        e = e.argument;
      if (e?.type !== 'CallExpression')
        return false;
      const name = calleeName(e);
      return VIEW_CALLS.has(name) || name === 'create' && e.callee.type === 'MemberExpression' &&
        nameOf(e.callee.object) === 'View';
    };

    /** 'view' (it answers a platform view), 'none' (nothing) or 'value' (a piece its caller places). */
    const answers = (fn) => {
      const declared = fn.returnType ? source.getText(fn.returnType) : '';
      if (VIEW_TYPE.test(declared))
        return 'view';
      const values = answered(fn);
      if (values.length === 0)
        return 'none';
      return values.some((v) => makesView(v) || v.type === 'Identifier' && makesView(valueOf(context, v)) ||
        v.type === 'MemberExpression' && nameOf(v) === 'view') ? 'view' : 'value';
    };

    const newOf = (cls) => cls.id ? calls.filter((c) => c.type === 'NewExpression' &&
      c.callee.type === 'Identifier' && c.callee.name === cls.id.name) : null;

    /** Where this file calls `fn`, or hands it on by name; null when nothing names it. */
    const sitesOf = (fn) => {
      const p = fn.parent;
      let plain = null;
      let member = null;
      if (fn.type === 'FunctionDeclaration' && fn.id)
        plain = fn.id.name;
      else if (p.type === 'VariableDeclarator' && p.init === fn && p.id.type === 'Identifier')
        plain = p.id.name;
      else if (p.type === 'MethodDefinition' && p.kind === 'constructor')
        return newOf(p.parent.parent);
      else if ((p.type === 'MethodDefinition' && p.kind === 'method' || p.type === 'PropertyDefinition' ||
        p.type === 'Property') && p.value === fn)
        member = keyName(p.key);
      if (plain === null && member === null)
        return null;
      const outside = (n) => n.range[0] < fn.range[0] || n.range[1] > fn.range[1];
      const sites = calls.filter((c) => outside(c) && (plain !== null ?
        c.callee.type === 'Identifier' && c.callee.name === plain :
        c.callee.type === 'MemberExpression' && !c.callee.computed && keyName(c.callee.property) === member));
      return plain === null ? sites : [...sites, ...handed.filter((r) => r.name === plain && outside(r))];
    };

    const callers = (sites, depth) => {
      if (sites === null || sites.length === 0)
        return 'unowned';
      if (depth >= DEPTH)
        return 'owned';
      return sites.every((s) => !['unowned', 'await'].includes(verdict(s, depth + 1))) ? 'owned' : 'unowned';
    };

    const inControlClass = (node) => {
      for (let p = node.parent; p; p = p.parent) {
        if (p.type === 'ClassBody' && controlClass(p.parent))
          return true;
      }
      return false;
    };

    const enclosingFunction = (node) => {
      let p = node.parent;
      while (p && !isFunction(p))
        p = p.parent;
      return p ?? null;
    };

    const functionName = (fn) => fn.type === 'FunctionDeclaration' ? fn.id?.name ?? null :
      fn.parent.type === 'VariableDeclarator' && fn.parent.init === fn && fn.parent.id.type === 'Identifier' ?
        fn.parent.id.name : null;

    /** What the name `id` holds where it is read, through its own scopes. */
    const bindingOf = (id) => unwrap(valueOf(context, id));

    /** `ui.dialog(…)`, a chain of its `.add(…)`s, or a name bound to one. */
    const platformDialog = (expr) => {
      let e = unwrap(expr);
      if (e?.type === 'Identifier')
        e = bindingOf(e);
      while (e?.type === 'CallExpression' && e.callee.type === 'MemberExpression' && calleeName(e) !== 'dialog')
        e = unwrap(e.callee.object);
      return e?.type === 'CallExpression' && calleeName(e) === 'dialog' && imports.resolve(e.callee) === null;
    };

    /** The message for a call that hands an element to the platform to show, keep and close; else null. */
    const hostKind = (call) => {
      const c = call.callee;
      if (c.type !== 'MemberExpression')
        return null;
      const name = calleeName(call);
      const on = nameOf(c.object);
      if (name === 'newView' || name === 'addView')
        return 'hostView';
      if (name === 'fromRoot')
        return on === 'View' ? 'hostView' : on === 'Widget' ? 'hostWidget' : null;
      if (name === 'dock')
        return on === 'dockManager' ? 'hostDock' : null;
      return name === 'add' && platformDialog(c.object) ? 'hostDialog' : null;
    };

    /** `<control>.root`, the control built here or by a name bound to one: the control and how to name it. */
    const rootOf = (expr) => {
      const e = unwrap(expr);
      if (e?.type !== 'MemberExpression' || e.computed || e.property.type !== 'Identifier' || e.property.name !== 'root')
        return null;
      const object = unwrap(e.object);
      if (built.includes(object))
        return {node: object, what: `${builtControl(object)}.root`};
      const bound = object.type === 'Identifier' ? bindingOf(object) : null;
      return bound !== null && built.includes(bound) ? {node: bound, what: `${object.name}.root`} : null;
    };

    /** 'owned', 'unowned', 'await' (built after an `await`), or 'class' (a control's own code: build-in-owner's). */
    function verdict(node, depth) {
      if (inControlClass(node))
        return 'class';
      let arg = node;
      while (WRAPPERS.has(arg.parent.type))
        arg = arg.parent;
      // handed to what takes it over: `host(control)`, `ViewHost.of(view).dock(control, …)`
      const taker = arg.parent;
      if (taker.type === 'CallExpression' && taker.arguments.includes(arg) &&
        (calleeName(taker) === 'dock' || calleeName(taker) !== 'add' && owns(taker)))
        return 'owned';
      for (let child = node, p = node.parent; p; child = p, p = p.parent) {
        if (child === node && node.type === 'Identifier') {
          if (owns(p))
            return 'owned';
          if (runsNow(p))
            continue;
          return 'unowned';
        }
        if (p.type === 'PropertyDefinition' && p.value === child && !isFunction(child))
          return p.static ? 'unowned' : callers(newOf(p.parent.parent), depth);
        if (!isFunction(p))
          continue;
        if (afterAwait(p, node))
          return 'await';
        const call = handedTo(p);
        if (call !== null && owns(call))
          return 'owned';
        if (returns(p, node))
          return 'owned';
        if (call !== null && runsNow(call)) {
          p = call;
          continue;
        }
        // a builder of a piece is its callers' to place and own; one that makes a view, or answers nothing, put
        // the control somewhere for good
        const sites = sitesOf(p);
        if ((sites === null || sites.length === 0) && answers(p) === 'value')
          return 'owned';
        return callers(sites, depth);
      }
      return 'unowned';
    }

    return {
      Program(node) {
        imports = new U2Imports(node);
      },
      CallExpression(node) {
        const c = node.callee;
        if (c.type === 'MemberExpression' && !c.computed && c.property.name === 'runWith' && isScope(c.object) &&
          isNothing(node.arguments[0]))
          context.report({node, messageId: 'runWith'});
        calls.push(node);
        handed.push(...node.arguments.filter((a) => a.type === 'Identifier'));
        if (builtControl(node) !== null)
          built.push(node);
      },
      NewExpression(node) {
        if (isScope(node.callee) && isNothing(node.arguments[0]))
          context.report({node, messageId: 'newScope'});
        calls.push(node);
        handed.push(...node.arguments.filter((a) => a.type === 'Identifier'));
        if (builtControl(node) !== null)
          built.push(node);
      },
      'Program:exit'() {
        for (let grew = true; grew;) {
          grew = false;
          for (const node of built) {
            const fn = enclosingFunction(node);
            const name = fn !== null && returns(fn, node) ? functionName(fn) : null;
            if (name !== null && !localBuilders.has(name)) {
              localBuilders.add(name);
              grew = true;
            }
          }
          if (grew)
            built.push(...calls.filter((c) => !built.includes(c) && builtControl(c) !== null));
        }
        const hosted = new Set();
        for (const call of calls) {
          const messageId = call.type === 'CallExpression' ? hostKind(call) : null;
          if (messageId === null)
            continue;
          for (const part of call.arguments.flatMap((a) => valueParts(a))) {
            const control = rootOf(part);
            if (control !== null) {
              hosted.add(control.node);
              context.report({node: part, messageId, data: {what: control.what}});
            }
          }
        }
        for (const node of built) {
          const found = hosted.has(node) ? 'hosted' : verdict(node, 0);
          if (found === 'unowned' || found === 'await') {
            context.report({node, messageId: found === 'await' ? 'afterAwait' :
              node.type === 'NewExpression' ? 'unownedControl' : 'unownedCall', data: {what: builtControl(node)}});
          }
        }
      },
    };
  },
};
