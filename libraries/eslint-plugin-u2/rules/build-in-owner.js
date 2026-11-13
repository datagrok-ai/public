'use strict';

const {DEFAULT_CONTROLS, calleeName, isFunction, sourceOf, walk} = require('../lib/u2-imports');

const nameOf = (node) => node.type === 'Identifier' ? node.name :
  node.type === 'MemberExpression' && !node.computed && node.property.type === 'Identifier' ? node.property.name : null;

const receiverOf = (call) => call.callee.type === 'MemberExpression' ? nameOf(call.callee.object) : null;

const ownedBuild = (call) => {
  const name = calleeName(call);
  const receiver = receiverOf(call);
  return name === 'runInScope' || name === 'build' && receiver !== 'Control' ||
    name === 'runWith' && receiver === 'Scope';
};

module.exports = {
  meta: {
    type: 'problem',
    docs: {description: 'A control builds its parts inside its own scope'},
    messages: {
      outside: '{{what}} in {{cls}}\'s {{where}} is built outside {{cls}}\'s scope, so it belongs to whoever is ' +
        'building {{cls}} and outlives it. Build it inside the control\'s scope: `this.build(() => …)`.',
    },
    schema: [{type: 'object', properties: {
      controls: {type: 'array', items: {type: 'string'}},
      factories: {type: 'array', items: {type: 'string'}},
      methods: {type: 'boolean'},
    }, additionalProperties: false}],
  },
  create(context) {
    const options = context.options[0] ?? {};
    const controls = options.controls ? new Set(options.controls) : null;
    const factories = new Set(options.factories ?? []);
    const isControl = (name) => name !== null && (controls === null ? DEFAULT_CONTROLS.test(name) : controls.has(name));
    const keys = sourceOf(context).visitorKeys;

    const check = (cls, part, where) => {
      walk(part, keys, (n) => {
        if (isFunction(n))
          return false;
        if (n.type === 'CallExpression' && ownedBuild(n))
          return false;
        const call = n.type === 'CallExpression';
        const name = n.type === 'NewExpression' || call ? nameOf(n.callee) : null;
        const controlBuild = call && name === 'build' && receiverOf(n) === 'Control';
        const built = n.type === 'NewExpression' ? isControl(name) : call && (controlBuild || factories.has(name));
        if (built) {
          const what = n.type === 'NewExpression' ? `new ${name}` : controlBuild ? 'Control.build' : name;
          context.report({node: n, messageId: 'outside', data: {what: `\`${what}(…)\``, cls, where}});
        }
        return !built;
      });
    };

    return {
      'ClassDeclaration, ClassExpression'(node) {
        if (!node.superClass || !isControl(nameOf(node.superClass)))
          return;
        const cls = node.id?.name ?? 'this control';
        for (const m of node.body.body) {
          if (m.type === 'MethodDefinition' && m.kind === 'constructor' && m.value.body)
            check(cls, m.value.body, 'constructor');
          else if (m.type === 'PropertyDefinition' && !m.static && m.value && !isFunction(m.value))
            check(cls, m.value, 'field initializer');
          else if (options.methods && m.type === 'MethodDefinition' && m.kind !== 'constructor' && !m.static &&
            m.value.body) {
            const kind = m.kind === 'method' ? 'method' : `${m.kind}ter`;
            check(cls, m.value.body, `${kind} ${nameOf(m.key) ?? ''}`.trim());
          }
        }
      },
    };
  },
};
