'use strict';

const ENTRIES = ['', '/dg', '/testing', '/filter', '/dg/filter', '/package.json'];

function entryFor(subpath) {
  if (/^\/src\/(core\/)?filter(\/|$)/.test(subpath))
    return '@datagrok-libraries/u2/filter';
  if (/^\/src\/dg\/filter(\/|$)/.test(subpath))
    return '@datagrok-libraries/u2/dg/filter';
  if (/^\/src\/dg(\/|$)/.test(subpath))
    return '@datagrok-libraries/u2/dg';
  if (/^\/tests(\/|$)/.test(subpath))
    return '@datagrok-libraries/u2/testing';
  return '@datagrok-libraries/u2';
}

module.exports = {
  meta: {
    type: 'problem',
    docs: {description: 'Import u2 through its entry points, never a path inside the package'},
    messages: {
      deep: '`{{source}}` reaches inside u2: only its entry points are a contract. Import from `{{entry}}`.',
    },
    schema: [{type: 'object', properties: {entries: {type: 'array', items: {type: 'string'}}},
      additionalProperties: false}],
  },
  create(context) {
    const entries = context.options[0]?.entries ?? ENTRIES;
    const check = (source) => {
      if (source?.type !== 'Literal' || typeof source.value !== 'string')
        return;
      const m = /^@datagrok-libraries\/u2(\/.*)?$/.exec(source.value);
      if (m === null)
        return;
      const subpath = m[1] ?? '';
      if (entries.includes(subpath) || /^\/css\/[^/]+\.css$/.test(subpath))
        return;
      context.report({node: source, messageId: 'deep', data: {source: source.value, entry: entryFor(subpath)}});
    };
    return {
      ImportDeclaration: (node) => check(node.source),
      ExportNamedDeclaration: (node) => check(node.source),
      ExportAllDeclaration: (node) => check(node.source),
      ImportExpression: (node) => check(node.source),
    };
  },
};
