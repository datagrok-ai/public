'use strict';

const {calleeName, isFunction, sourceOf, staticText, walk} = require('../lib/u2-imports');

const EMPTY = /^(disabled|not available|unavailable|n\/a|not applicable)\.?$/i;
const REASON_KEYS = new Set(['disabled', 'disabledReason', 'reason']);
const REASON_CALLS = new Set(['okDisabled', 'setDisabled', 'disable']);
const WRAPPERS = new Set(['computed', 'signal']);

const keyName = (node) => node.type === 'Identifier' ? node.name :
  node.type === 'Literal' && typeof node.value === 'string' ? node.value : null;

module.exports = {
  meta: {
    type: 'suggestion',
    docs: {description: 'A control that is off says why, in words a user can act on'},
    messages: {
      empty: '"{{text}}" says the control is off, not why. Name what turns it on: "Pick a source file first", ' +
        '"Only the campaign\'s owner can edit it".',
    },
    schema: [],
  },
  create(context) {
    const keys = sourceOf(context).visitorKeys;
    const results = (node, out = []) => {
      if (node.type === 'ConditionalExpression') {
        results(node.consequent, out);
        results(node.alternate, out);
      }
      else if (node.type === 'LogicalExpression') {
        results(node.left, out);
        results(node.right, out);
      }
      else if (node.type === 'TSAsExpression' || node.type === 'TSNonNullExpression')
        results(node.expression, out);
      else if (isFunction(node)) {
        if (node.body.type !== 'BlockStatement')
          results(node.body, out);
        else {
          walk(node.body, keys, (n) => {
            if (isFunction(n))
              return false;
            if (n.type === 'ReturnStatement' && n.argument)
              results(n.argument, out);
            return true;
          });
        }
      }
      else if (node.type === 'CallExpression' && WRAPPERS.has(calleeName(node)) && node.arguments[0])
        results(node.arguments[0], out);
      else if (node.type === 'Literal' || node.type === 'TemplateLiteral' && node.expressions.length === 0)
        out.push(node);
      return out;
    };
    const check = (value) => {
      for (const node of results(value)) {
        const text = staticText(node)?.trim();
        if (text !== undefined && EMPTY.test(text))
          context.report({node, messageId: 'empty', data: {text}});
      }
    };
    return {
      Property(node) {
        if (!node.computed && REASON_KEYS.has(keyName(node.key)) && node.parent.type === 'ObjectExpression')
          check(node.value);
      },
      AssignmentExpression(node) {
        const target = node.left;
        if (target.type === 'MemberExpression' && !target.computed && REASON_KEYS.has(target.property.name))
          check(node.right);
      },
      CallExpression(node) {
        const name = calleeName(node);
        if (REASON_CALLS.has(name) && node.arguments[0])
          check(node.arguments[0]);
        else if (name === 'setAttribute' && staticText(node.arguments[0] ?? {type: 'none'}) === 'aria-description' &&
          node.arguments[1])
          check(node.arguments[1]);
      },
    };
  },
};
