'use strict';

const {INTERNAL_CLASS, INTERNAL_SELECTOR, calleeName, staticText} = require('../lib/u2-imports');

const SELECTOR_CALLS = new Set(['querySelector', 'querySelectorAll', 'closest', 'matches']);
const CLASS_CALLS = new Set(['add', 'remove', 'toggle', 'contains', 'replace']);

const onClassList = (call) => call.callee.type === 'MemberExpression' &&
  call.callee.object.type === 'MemberExpression' && !call.callee.object.computed &&
  call.callee.object.property.name === 'classList';

module.exports = {
  meta: {
    type: 'problem',
    docs: {description: 'No u2 class names or data-u2 attributes in consumer lookups'},
    messages: {
      internal: '`{{text}}` reaches into u2 internals; its classes and data-u2 attributes are not a contract. Keep ' +
        'a reference to the control and read its state, or find it by automation name with the u2 testing kit.',
    },
    schema: [],
  },
  create(context) {
    return {
      CallExpression(node) {
        const name = calleeName(node);
        if (node.callee.type !== 'MemberExpression' || name === null)
          return;
        if (SELECTOR_CALLS.has(name)) {
          const text = node.arguments[0] && staticText(node.arguments[0]);
          if (text && INTERNAL_SELECTOR.test(text))
            context.report({node: node.arguments[0], messageId: 'internal', data: {text}});
        }
        else if (name === 'getElementsByClassName' || CLASS_CALLS.has(name) && onClassList(node)) {
          for (const arg of node.arguments) {
            const text = staticText(arg);
            if (text && INTERNAL_CLASS.test(text))
              context.report({node: arg, messageId: 'internal', data: {text}});
          }
        }
      },
    };
  },
};
