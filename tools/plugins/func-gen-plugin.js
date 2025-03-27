const fs = require('fs');
const path = require('path');
const {parse} = require('@babel/parser');
const traverse = require('@babel/traverse').default;
const {reservedDecorators, getFuncAnnotation, generateImport, generateExport} = require('../bin/utils/func-generation');


class FuncGeneratorPlugin {
  constructor(options = {outputPath: './src/package.g.ts'}) {
    this.options = options;
  }

  apply(compiler) {
    const srcDirPath = path.join(compiler.context, 'src');
    let packageFilePath = path.join(srcDirPath, 'package.ts');
    if(!fs.existsSync(packageFilePath))
      packageFilePath = path.join(srcDirPath, 'package.js');
    const tsFiles = this._getTsFiles(srcDirPath);
    const genImports = [];
    const genExports = [];

    compiler.hooks.compilation.tap('FuncGeneratorPlugin', (_compilation) => {
      this._clearGeneratedFile();

      for (const file of tsFiles) {
        const content = fs.readFileSync(file, 'utf-8');
        if (!content)
          continue;
        const ast = parse(content, {
          sourceType: 'module',
          plugins: ['typescript', [
            'decorators',
            {decoratorsBeforeExport: true},
          ]],
        });
        const functions = [];
        const imports = [];
  
        traverse(ast, {
          ClassDeclaration: (nodePath) => this._addNodeData(nodePath.node, file,
            srcDirPath, functions, imports, genImports, genExports),
        });
  
        if (fs.existsSync(this.options.outputPath)) { 
          const content = fs.readFileSync(this.options.outputPath, 'utf-8');
          const output = content ? this._insertImports(content, imports) : imports.join('\n');
          fs.writeFileSync(this.options.outputPath, output, 'utf-8');
        } else 
          fs.writeFileSync(this.options.outputPath, imports.join('\n'), 'utf-8');
        
  
        fs.appendFileSync(this.options.outputPath, functions.join('\n'), 'utf-8');
      }

      this._writeToPackageFile(packageFilePath, genImports, genExports);
    });
  }

  _getTsFiles(root) {
    const tsFiles = [];
    const extPattern = /\.tsx?$/;
    const excludedFiles = ['package.ts', 'package-test.ts', 'package.g.ts'];

    function findFiles(dir) {
      const files = fs.readdirSync(dir);
      for (const file of files) {
        const fullPath = path.join(dir, file);
        if (fs.statSync(fullPath).isDirectory()) 
          findFiles(fullPath);
        else if (extPattern.test(file) && !excludedFiles.includes(file)) 
          tsFiles.push(fullPath);
        
      }
    }

    findFiles(root);
    return tsFiles;
  }

  _clearGeneratedFile() {
    fs.writeFileSync(this.options.outputPath, '');
  }

  _writeToPackageFile(filePath, imports, exports) {
    if (imports.length !== exports.length)
      return;
    let content = fs.readFileSync(filePath, 'utf-8');
    for (let i = 0; i < imports.length; i++) {
      const importStatement = imports[i];
      const exportStatement = exports[i];
      if (!content.includes(importStatement.trim()))
        content = importStatement + content + exportStatement;
    }
    fs.writeFileSync(filePath, content, 'utf-8');
  }

  _addNodeData(node, file, srcDirPath, functions, imports, genImports, genExports) {
    if (!node.decorators || !node.decorators.length)
      return;

    function modifyImportPath(dirPath, filePath) {
      const relativePath = path.relative(dirPath, filePath);
      return `./${relativePath.slice(0, relativePath.length - 3).replace(/\\/g, '/')}`;
    }

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
        imports.push(generateImport(className, modifyImportPath(path.dirname(this.options.outputPath), file)));
        const funcName = `_${className}`;
        functions.push(reservedDecorators[name]['genFunc'](getFuncAnnotation({
          ...reservedDecorators[name]['metadata'],
          ...options,
        }), className));

        genImports.push(generateImport(funcName, modifyImportPath(srcDirPath, this.options.outputPath)));
        genExports.push(generateExport(funcName));
      }
    }
  }

  _insertImports(content, imports) {
    const ast = parse(content, {sourceType: 'module', plugins: ['typescript', 'decorators']});
    let lastImportLoc = null;
    traverse(ast, {
      ImportDeclaration(nodePath) {
        lastImportLoc = nodePath.node.end + 1;
      },
    });
    return lastImportLoc === null ? imports.join('\n') + content :
      content.slice(0, lastImportLoc) + imports.join('\n') +
      (content[lastImportLoc] === '\n' ? '' : '\n') +
      content.slice(lastImportLoc, content.length);
  }
}

module.exports = FuncGeneratorPlugin;
