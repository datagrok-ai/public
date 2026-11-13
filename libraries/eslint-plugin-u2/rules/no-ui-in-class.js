'use strict';

const {U2Imports, calleeName, isFunction, sourceOf, walk} = require('../lib/u2-imports');

const RUNS_NOW = new Set(['runInScope', 'build', 'runWith']);

function constructionParts(cls) {
  const parts = [];
  for (const m of cls.body.body) {
    if (m.type === 'MethodDefinition' && m.kind === 'constructor' && m.value.body)
      parts.push(m.value.body);
    else if (m.type === 'PropertyDefinition' && !m.static && m.value && !isFunction(m.value))
      parts.push(m.value);
  }
  return parts;
}

function rootOf(node) {
  let n = node;
  while (n.type === 'MemberExpression' || n.type === 'CallExpression')
    n = n.type === 'MemberExpression' ? n.object : n.callee;
  return n;
}

const runsNow = (fn) => fn.parent.type === 'CallExpression' && fn.parent.arguments.includes(fn) &&
  RUNS_NOW.has(calleeName(fn.parent));

module.exports = {
  meta: {
    type: 'problem',
    docs: {description: 'A consumer class does not build u2 UI while it is constructed'},
    messages: {
      ui: '{{name}} builds u2 UI ({{what}}) while it is constructed, so the UI belongs to whoever is building it and ' +
        'outlives the object. Make it a function component: `component((self) => ({root, ...api}), {name})`.',
    },
    schema: [{type: 'object', properties: {
      controls: {type: 'array', items: {type: 'string'}},
      factories: {type: 'array', items: {type: 'string'}},
    }, additionalProperties: false}],
  },
  create(context) {
    const options = context.options[0] ?? {};
    const controls = options.controls ? new Set(options.controls) : null;
    const factories = new Set(options.factories ?? []);
    const keys = sourceOf(context).visitorKeys;
    let imports;

    const isControl = (name) => controls === null ? /^[A-Z]/.test(name) : controls.has(name);
    const ownsItself = (call) => {
      const root = imports.resolve(rootOf(call.callee));
      return root !== null && (root.name === 'ViewHost' || root.name === 'Control' && calleeName(call) === 'build');
    };
    const builtUi = (part) => {
      let what = null;
      walk(part, keys, (n) => {
        if (what !== null)
          return false;
        if (isFunction(n))
          return runsNow(n);
        if (n.type === 'CallExpression' && ownsItself(n))
          return false;
        const u2 = n.type === 'NewExpression' || n.type === 'CallExpression' ? imports.resolve(n.callee) : null;
        if (u2 !== null && n.type === 'NewExpression' && isControl(u2.name))
          what = `new ${u2.name}(…)`;
        else if (u2 !== null && n.type === 'CallExpression' && factories.has(u2.name))
          what = `${u2.name}(…)`;
        return true;
      });
      return what;
    };

    return {
      Program(node) {
        imports = new U2Imports(node);
      },
      'ClassDeclaration, ClassExpression'(node) {
        if (node.superClass && imports.resolve(node.superClass))
          return;
        for (const part of constructionParts(node)) {
          const what = builtUi(part);
          if (what !== null) {
            context.report({node: node.id ?? node, messageId: 'ui', data: {name: node.id?.name ?? 'This class', what}});
            return;
          }
        }
      },
    };
  },
};
