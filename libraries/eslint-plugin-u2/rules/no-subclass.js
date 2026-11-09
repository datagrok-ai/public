'use strict';

const {U2Imports} = require('../lib/u2-imports');

module.exports = {
  meta: {
    type: 'problem',
    docs: {description: 'Consumer code writes function components; a subclass of a u2 or ems class records why'},
    messages: {
      subclass: '{{name}} extends {{lib}}\'s {{base}}. Function components are the documented form: return ' +
        '`component((self) => ({root, ...api}), {name})` from a function, and keep statics on a holder class that ' +
        'builds no UI. A class that is truly needed records its reason: ' +
        '`node libraries/u2/scripts/check-consumer.mjs --update-baseline --allow-growth "<why>"`.',
    },
    schema: [],
  },
  create(context) {
    let imports;
    return {
      Program(node) {
        imports = new U2Imports(node);
      },
      'ClassDeclaration, ClassExpression'(node) {
        const base = node.superClass && imports.resolve(node.superClass);
        if (base)
          context.report({node: node.id ?? node, messageId: 'subclass',
            data: {name: node.id?.name ?? 'This class', lib: base.lib, base: base.name}});
      },
    };
  },
};
