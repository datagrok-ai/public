'use strict';

const {U2Imports, isFunction, sourceOf} = require('../lib/u2-imports');

const isHandler = (node) => isFunction(node) || node.type === 'ObjectExpression';

module.exports = {
  meta: {
    type: 'problem',
    docs: {description: 'In code that builds u2 UI, every subscription has an owner'},
    messages: {
      raw: 'A raw `.subscribe(…)` is owned by nothing and outlives the control unless unsubscribed by hand. ' +
        'Subscribe on the owner: `self.subscribe(source, handler)` in a component, `host.subscribe(source, handler)` ' +
        'on a view\'s `ViewHost`, `owner.subscribe(source, handler)` on any other owner. It ends with the owner and ' +
        'contains the handler.',
    },
    schema: [],
  },
  create(context) {
    if (!U2Imports.usesU2(sourceOf(context).text))
      return {};
    return {
      CallExpression(node) {
        const c = node.callee;
        if (c.type !== 'MemberExpression' || c.computed || c.property.name !== 'subscribe')
          return;
        if (node.arguments.length < 2 || isHandler(node.arguments[0]))
          context.report({node: c.property, messageId: 'raw'});
      },
    };
  },
};
