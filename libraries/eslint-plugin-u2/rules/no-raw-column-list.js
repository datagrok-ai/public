'use strict';

const {U2Imports, sourceOf, walk} = require('../lib/u2-imports');

const CHOICE_KEYS = new Set(['items', 'choices']);
const LISTINGS = new Set(['toList', 'names', 'byNames', 'byTags', 'categorical', 'numerical', 'dateTime', 'boolean']);

const keyName = (node) => node.type === 'Identifier' ? node.name :
  node.type === 'Literal' && typeof node.value === 'string' ? node.value : null;

/** `x.columns` read as a list: `.toList()`, `.names()` and the like, `Array.from(x.columns)`, `[...x.columns]`. */
function listing(columns) {
  const parent = columns.parent;
  if (parent.type === 'MemberExpression' && parent.object === columns && !parent.computed &&
    LISTINGS.has(parent.property.name))
    return `.columns.${parent.property.name}${parent.parent.type === 'CallExpression' ? '()' : ''}`;
  if (parent.type === 'SpreadElement')
    return '[...columns]';
  if (parent.type !== 'CallExpression' || parent.arguments[0] !== columns)
    return null;
  const callee = parent.callee;
  if (callee.type === 'Identifier' && callee.name === 'wu')
    return 'wu(….columns)';
  return callee.type === 'MemberExpression' && callee.object.name === 'Array' && callee.property.name === 'from' ?
    'Array.from(….columns)' : null;
}

module.exports = {
  meta: {
    type: 'problem',
    docs: {description: 'A u2 choice of columns lists them through the platform, not from `df.columns`'},
    messages: {
      raw: '`{{key}}` gets a column list built here ({{what}}), which offers the hidden `~` columns and nameless ' +
        'ones. Offer columns with `columnInput(label, table, {filter})`, or list them with ' +
        '`Platform.listedColumns(table, filter)`, so the platform\'s rule applies.',
    },
    schema: [],
  },
  create(context) {
    const source = sourceOf(context);
    const keys = source.visitorKeys;
    let imports;

    const ofU2Call = (property) => {
      for (let at = property.parent; at; at = at.parent) {
        if (at.type === 'CallExpression' || at.type === 'NewExpression')
          return imports.resolve(at.callee) !== null;
        if (!['ObjectExpression', 'ArrayExpression', 'Property', 'SpreadElement'].includes(at.type))
          return false;
      }
      return false;
    };

    const rawIn = (value, seen) => {
      let found = null;
      walk(value, keys, (n) => {
        if (found !== null)
          return false;
        if (n.type === 'MemberExpression' && !n.computed && n.property.name === 'columns')
          found = listing(n);
        else if (n.type === 'Identifier' && !seen.has(n.name)) {
          const variable = source.getScope(n).references.find((r) => r.identifier === n)?.resolved;
          const init = variable?.defs[0]?.node.type === 'VariableDeclarator' ? variable.defs[0].node.init : null;
          if (init)
            found = rawIn(init, new Set([...seen, n.name]));
        }
        return found === null;
      });
      return found;
    };

    return {
      Program(node) {
        imports = new U2Imports(node);
      },
      Property(node) {
        const key = node.computed ? null : keyName(node.key);
        if (!CHOICE_KEYS.has(key) || node.parent.type !== 'ObjectExpression' || !ofU2Call(node))
          return;
        const what = rawIn(node.value, new Set());
        if (what !== null)
          context.report({node: node.value, messageId: 'raw', data: {key, what}});
      },
    };
  },
};
