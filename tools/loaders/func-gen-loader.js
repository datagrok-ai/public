const { parse } = require('@babel/parser');
const traverse = require('@babel/traverse').default;
const fs = require('fs');
const path = require('path');

// Works if there's a generated .js file (run `tsc constants.ts`)
const reservedDecorators = require('./constants').reservedDecorators;
const getFuncAnnotation = require('./constants').getFuncAnnotation;
const generateImport = require('./constants').generateImport;


module.exports = function loader(source) {
  const ast = parse(source, { sourceType: 'module', plugins: ['typescript', ['decorators', {decoratorsBeforeExport: true}]] });
  const functions = [];
  const imports = [];
  const resourcePath = this.resourcePath;
  const outputPath = 'package.g.ts';

  traverse(ast, {
    ClassDeclaration(nodePath) {
      const { node } = nodePath;
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

            const className = node.id.name;
            const filePath = path.relative(process.cwd(), resourcePath);
            const modifiedPath = `./${filePath.slice(0, filePath.length - 3).replace(/\\/g, '/')}`;
            imports.push(generateImport(className, modifiedPath));
            functions.push(reservedDecorators[name]['genFunc'](getFuncAnnotation({
              ...reservedDecorators[name]['metadata'],
              ...options,
            }), className));
          }
        }
      }
    }
  });

  if (fs.existsSync(outputPath)) {
    function insertImports(content) {
      const ast = parse(content, { sourceType: 'module', plugins: ['typescript', 'decorators'] });
      let lastImportLoc = null;
      traverse(ast, {
        ImportDeclaration(nodePath) {
          lastImportLoc = nodePath.node.end + 1;
        }
      });
      return lastImportLoc === null ? imports.join('\n') + content
        : content.slice(0, lastImportLoc) + imports.join('\n') +
        (content[lastImportLoc] === '\n' ? '' : '\n') +
        content.slice(lastImportLoc, content.length);
    }

    const content = fs.readFileSync(outputPath, 'utf-8');
    const output = content ? insertImports(content) : imports.join('\n');
    fs.writeFileSync(outputPath, output, 'utf-8');
  } else {
    fs.writeFileSync(outputPath, imports.join('\n'), 'utf-8');
  }

  fs.appendFileSync(outputPath, functions.join('\n'), 'utf-8');
  return source;
}
