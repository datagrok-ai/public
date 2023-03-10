const { parse } = require('@babel/parser');
const traverse = require('@babel/traverse').default;
const fs = require('fs');

// Works if there's a generated .js file (run `tsc constants.ts`)
const reservedDecorators = require('./constants').reservedDecorators;
const getFuncAnnotation = require('./constants').getFuncAnnotation;


module.exports = function loader(source) {
  const ast = parse(source, { sourceType: 'module', plugins: ['typescript', 'decorators'] });
  const functions = [];

  traverse(ast, {
    ClassDeclaration(path) {
      const { node } = path;
      if (node.decorators && node.decorators.length > 0) {
        for (const decorator of node.decorators) {
          const exp = decorator.expression;
          const name = exp.callee.property.name;
          const options = {};
          if (name in reservedDecorators) {
            if (exp.arguments && exp.arguments.length === 1) {
              const props = exp.arguments[0].properties;
              for (const prop of props)
                options[prop.key.name] = prop.value.value;
            }

            functions.push(reservedDecorators[name]['genFunc'](getFuncAnnotation({
              ...reservedDecorators[name]['metadata'],
              ...options,
            }), node.id.name));
          }
        }
      }
    }
  });

  fs.appendFileSync('package.g.ts', functions.join('\n'), 'utf-8');

  return JSON.stringify(source);
}
