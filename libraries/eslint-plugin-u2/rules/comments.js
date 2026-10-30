'use strict';

const {sourceOf} = require('../lib/u2-imports');

const DIRECTIVE = /^[\s*]*(eslint|global|exported|istanbul|c8|@ts-|prettier|webpack|jshint|@license|copyright\b)/i;
const MEMBERS = ['PropertyDefinition', 'MethodDefinition', 'AccessorProperty', 'TSAbstractPropertyDefinition',
  'TSAbstractMethodDefinition', 'TSAbstractAccessorProperty', 'TSParameterProperty'];

function access(member) {
  if (member.accessibility === 'private' || member.accessibility === 'protected')
    return member.accessibility;
  return member.key?.type === 'PrivateIdentifier' ? 'private' : null;
}

module.exports = {
  meta: {
    type: 'suggestion',
    docs: {description: 'No file-header narratives, and doc blocks only on the API a caller reads'},
    messages: {
      header: 'A file header narrates what the code, the README or the history already say, and goes stale. ' +
        'Delete it; a one-line `//` why stays next to the code it explains.',
      privateDoc: 'A doc block on a {{access}} member documents no API. Delete it, or keep a one-line `//` why ' +
        'where the code cannot say it.',
    },
    schema: [],
  },
  create(context) {
    const source = sourceOf(context);
    return {
      Program(program) {
        const first = source.getAllComments()[0];
        const code = source.getFirstToken(program);
        if (first?.type === 'Block' && !first.value.startsWith('*') && !first.value.startsWith('!') &&
          !DIRECTIVE.test(first.value) && (code == null || first.range[1] <= code.range[0]))
          context.report({loc: first.loc, messageId: 'header'});
      },
      [MEMBERS.join(', ')](member) {
        const level = access(member);
        if (level === null)
          return;
        const doc = source.getCommentsBefore(member).pop();
        if (doc?.type === 'Block' && doc.value.startsWith('*'))
          context.report({loc: doc.loc, messageId: 'privateDoc', data: {access: level}});
      },
    };
  },
};
