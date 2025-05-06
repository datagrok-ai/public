const fs = require("fs");
const path = require("path");

const tsParser = require("@typescript-eslint/typescript-estree");
const generate = require('@babel/generator').default;

const {
  reservedDecorators,
  getFuncAnnotation,
  generateImport,
  generateExport,
  typesToAnnotation
} = require("../bin/utils/func-generation");

const baseImport = '\nimport * as DG from \'datagrok-api/dg\';\n\n';

class FuncGeneratorPlugin {
  constructor(options = { outputPath: "./src/package.g.ts" }) {
    this.options = options;
  }

  apply(compiler) {
    const srcDirPath = path.join(compiler.context, "src");
    let packageFilePath = path.join(srcDirPath, "package.ts");
    if (!fs.existsSync(packageFilePath))
      packageFilePath = path.join(srcDirPath, "package.js");
    const tsFiles = this._getTsFiles(srcDirPath);
    const genImports = [];
    const genExports = [];

    compiler.hooks.compilation.tap("FuncGeneratorPlugin", (_compilation) => {
      this._clearGeneratedFile();

      for (const file of tsFiles) {
        const content = fs.readFileSync(file, "utf-8");
        if (!content) continue;
        const ast = tsParser.parse(content, {
          sourceType: "module",
          plugins: [
            ["decorators", { decoratorsBeforeExport: true }],
            "classProperties",
            "typescript",
          ],
        });
        const functions = [];
        let imports = new Set();

        const fileName = path.basename(file);
        if(fileName === 'package-functions.ts')
          imports = new Set([... this._readImports(content)]);

        this._walk(ast, (node, parentClass) => {
          const decorators = node.decorators;
          if (!decorators || decorators.length === 0) return;

          if (node?.type === "ClassDeclaration") 
            this._addNodeData(node, file, srcDirPath, functions, imports, genImports, genExports);

          if (node?.type === "MethodDefinition")
            this._addNodeData(node, file, srcDirPath, functions, imports, genImports, genExports, parentClass);
          });
        this._insertImports([...imports]);
        fs.appendFileSync(this.options.outputPath, functions.join("\n"), "utf-8");
      }

      this._writeToPackageFile(packageFilePath, genImports, genExports);
    });
  }

  _addNodeData(node, file, srcDirPath, functions, imports, genImports, genExports, parent = undefined) {
    if (!node.decorators || !node.decorators.length || node.decorators?.length === 0)
      return;

    function modifyImportPath(dirPath, filePath) {
      const relativePath = path.relative(dirPath, filePath);
      return `./${relativePath.slice(0, relativePath.length - 3).replace(/\\/g, '/')}`;
    }

    const decorator =  node.decorators[0];
    const exp  = decorator.expression;
    const name = exp.callee?.property?.name || exp.callee?.name;
    const identifierName = node.id?.name || node.key?.name;
    const className = parent?.id?.name || parent?.key?.name;

    if (!name) return;

    const decoratorOptions = this._readDecoratorOptions(exp.arguments[0].properties);
    decoratorOptions.set('tags', [...(reservedDecorators[name]['metadata']['tags'] ?? [] ), ...(decoratorOptions.get('tags') ?? [])])
    const functionParams = node?.type === 'MethodDefinition' ? this._readMethodParamas(node) : [];

    const annotationByReturnType = node?.type === 'MethodDefinition' ? this._readFunctionReturnTypeInfo(node) : '';
    const annotationByReturnTypeObj = { name: 'result', type: annotationByReturnType };
    
    let importString = generateImport(node?.type === 'MethodDefinition' ? className : identifierName, modifyImportPath(path.dirname(this.options.outputPath), file));
    imports.add(importString);
    const funcName = `_${identifierName}`;
    const funcAnnotaionOptions = {
      ...reservedDecorators[name]['metadata'],
      ...(annotationByReturnType ? {outputs: [annotationByReturnTypeObj ?? {}]} : {}),
      ...Object.fromEntries(decoratorOptions),
      ...{inputs: functionParams},
    };    

    functions.push(reservedDecorators[name]['genFunc'](getFuncAnnotation(funcAnnotaionOptions), identifierName,'\n', (className ?? ''), functionParams));

    genImports.push(generateImport(funcName, modifyImportPath(srcDirPath, this.options.outputPath)));
    genExports.push(generateExport(funcName));
  }

  _readImports (content) {
    const importRegex = /(import(?:[\s\w{},*]+from\s*)?['"]([^'"]+)['"];)/g;
    const results = [];
  
    let match;
    while ((match = importRegex.exec(content)) !== null) {
      results.push(`${match[1]}\n`);
    }
    return results;
  }

  _readDecoratorOptions(properties){
    const resultMap = new Map();

    for(let prop of properties)
      resultMap.set(prop.key.name, this._evalLiteral(prop.value));

    return resultMap;
  }

  _evalLiteral(node) {
    if (!node) return null;

    switch (node.type) {
      case 'Literal':
        return node.value;

      case 'ArrayExpression':
        return node.elements.map(el => this._evalLiteral(el));

      case 'ObjectExpression':
        return Object.fromEntries(
          node.properties.map(p => {
            return [p.key.name || p.key.value, this._evalLiteral(p.value)];
          })
        );

      default:
        return '';
    }
  }

  _readMethodParamas(node) {
    const params = node?.value?.params?.map(param => {
      let baseParam = param.type === 'TSNonNullExpression' ? param.expression : param;
      let defaultValue = undefined;
      const options   = param.decorators?.length > 0?  Object.fromEntries(this._readDecoratorOptions(param.decorators[0]?.expression?.arguments[0].properties)) : undefined;
      
      if(baseParam.type === 'AssignmentPattern')
      {
        if (baseParam?.right?.type  === 'Literal')
          defaultValue = baseParam?.right?.raw;
        baseParam = baseParam?.left;
      } 

      if (baseParam.type === 'RestElement' || baseParam.type === 'Identifier') {
        let name =
          baseParam.type === 'RestElement'
            ? `...${baseParam.argument.name}`
            : baseParam.name;
    
        if (baseParam?.argument?.typeAnnotation)
          name += ': ' + generate(baseParam.argument.typeAnnotation.typeAnnotation).code;
    
        let type = '';
        if (baseParam?.typeAnnotation)
          type = generate(baseParam.typeAnnotation.typeAnnotation).code;
        else
          type = 'any';
    
        return { name: name, type: type, defaultValue: defaultValue, options: options };
      }
      else if (baseParam.type === 'ObjectPattern' || baseParam.type === 'ArrayPattern') {
        let name = '';
        if (baseParam.type === 'ObjectPattern') {
          const properties = baseParam.properties.map(prop => {
            if (prop.type === 'Property' && prop.key.type === 'Identifier')
              return prop.key.name;
            else if (prop.type === 'RestElement' && prop.argument.type === 'Identifier')
              return `...${prop.argument.name}`;
            else
              return generate(prop).code;
          });
          name = `{ ${properties.join(', ')} }`;
        } else {
          const elements = baseParam.elements.map(elem => {
            if (elem) {
              if (elem.type === 'Identifier')
                return elem.name;
              else
                return generate(elem).code;
            } else return '';
          });
          name = `[${elements.join(', ')}]`;
        }
    
        let type = '';
        if (baseParam.typeAnnotation)
          type = generate(baseParam.typeAnnotation.typeAnnotation).code;
        else
          type = 'any';

        return { name: name, type: type, defaultValue: defaultValue, options: options };
      }
      return { name: 'value', type: 'any', options: undefined };
    });
    return params;
  }

  _getTsFiles(root) {
    const tsFiles = [];
    const extPattern = /\.tsx?$/;
    const excludedFiles = ["package.ts", "package-test.ts", "package.g.ts"];

    function findFiles(dir) {
      const files = fs.readdirSync(dir);
      for (const file of files) {
        const fullPath = path.join(dir, file);
        if (fs.statSync(fullPath).isDirectory()) findFiles(fullPath);
        else if (extPattern.test(file) && !excludedFiles.includes(file))
          tsFiles.push(fullPath);
      }
    }

    findFiles(root);
    return tsFiles;
  }

  _walk(node, visitor, parent = null) {
    if (!node || typeof node !== "object") return;
  
    visitor(node, parent);
  
    for (const key in node) {
      const value = node[key];
      if (Array.isArray(value)) {
        value.forEach((child) => {
          if (child && typeof child.type === "string") {
            this._walk(child, visitor, node.type === 'ClassDeclaration'? node : parent );
          }
        });
      } else if (value && typeof value.type === "string") {
        this._walk(value, visitor, node.type === 'ClassDeclaration'? node : parent );
      }
    }
  }

  _readFunctionReturnTypeInfo(node) { 
    let resultType = 'any';
    let isArray = false;
    const annotation = node.value.returnType.typeAnnotation;
    if (annotation?.typeName?.name === 'Promise')
    {
      const argumnets = annotation.typeArguments?.params;
      if (argumnets && argumnets.length===1)
      {
        if (argumnets[0].typeName)
          resultType = argumnets[0].typeName?.right?.name ?? argumnets[0].typeName?.name;
        else if (argumnets[0].type !== 'TSArrayType')
          resultType = this._getTypeNameFromNode(argumnets[0]);
        else if (argumnets[0].elementType.type !== 'TSTypeReference'){
          isArray = true;
          resultType = this._getTypeNameFromNode(argumnets[0]?.elementType);
        }
        else{
          isArray = true;
          resultType = argumnets[0].elementType?.typeName?.name || argumnets[0].elementType?.typeName?.right?.name;
        }
      }
    }    
    else{      
      if (annotation.type === 'TSTypeReference')
        resultType = annotation.typeName?.right?.name ?? annotation.typeName?.name;
      else if (annotation.type !== 'TSArrayType')
        resultType = this._getTypeNameFromNode(annotation);
      else if  (annotation.elementType.type !== 'TSTypeReference'){
        isArray = true;
        resultType = this._getTypeNameFromNode(annotation?.elementType);
      }
      else{
        isArray = true;
        resultType = (annotation?.elementType?.typeName?.name || annotation?.elementType?.typeName?.right?.name);
      }
    }
    
    resultType = typesToAnnotation[resultType];
    if (isArray && resultType)
      resultType = `list<${resultType}>`
    return resultType;
  }

  _getTypeNameFromNode(typeNode) {
    if (typeNode.type === 'TSTypeReference') {
      return typeNode.typeName.name;
    } else if (typeNode.type === 'TSVoidKeyword') {
      return 'void';
    } else if (typeNode.type === 'TSNumberKeyword') {
      return 'number';
    } else if (typeNode.type === 'TSStringKeyword') {
      return 'string';
    } else if (typeNode.type === 'TSBooleanKeyword') {
      return 'boolean';
    } else {
      return typeNode.type;
    }
  }

  _clearGeneratedFile() {
    fs.writeFileSync(this.options.outputPath, baseImport);
  }

  _writeToPackageFile(filePath, imports, exports) {
    if (imports.length !== exports.length) return;
    let content = fs.readFileSync(filePath, "utf-8");
    for (let i = 0; i < imports.length; i++) {
      const importStatement = imports[i];
      const exportStatement = exports[i];
      if (!content.includes(importStatement.trim()))
        content = importStatement + content + exportStatement;
    }
    fs.writeFileSync(filePath, content, "utf-8");
  }

  _insertImports(importArray) {
    if (fs.existsSync(this.options.outputPath)) {
      const content = fs.readFileSync(this.options.outputPath, "utf-8");
      if (content)
        importArray.push(content);
      const output = importArray.join("");
      fs.writeFileSync(this.options.outputPath, `${output}`, "utf-8");
    } 
    else
      fs.writeFileSync(this.options.outputPath, `${baseImport}\n${importArray.join("")}`, "utf-8");
  }
}

module.exports = FuncGeneratorPlugin;