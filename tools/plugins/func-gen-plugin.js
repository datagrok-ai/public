const fs = require('fs');
const path = require('path');

const tsParser = require('@typescript-eslint/typescript-estree');
const generate = require('@babel/generator').default;

const {
  reservedDecorators,
  getFuncAnnotation,
  generateImport,
  generateExport,
  typesToAnnotation,
  dgAnnotationTypes,
  decoratorOptionToAnnotation,
  inputOptionsNames,
} = require('../bin/utils/func-generation');

const {toCamelCase} = require('../bin/commands/migrate');

const {api} = require('../bin/commands/api');

// Prebuilt CJS bundle of diff-grok's IVP parser (`getIVP` + a couple constants), tree-shaken
// to exclude the script-code generator. Regenerate with `npm run update:ivp-parser`.
function loadIvpParser() {
  try {
    const parser = require('./ivp-parser.bundle.cjs');
    if (typeof parser.getIVP !== 'function') {
      console.warn('[func-gen] ivp-parser.bundle.cjs has no getIVP export — skipping ivp models');
      return null;
    }
    return parser;
  } catch (e) {
    console.warn(`[func-gen] could not load ivp-parser.bundle.cjs (run "npm run update:ivp-parser"): ${e.message}`);
    return null;
  }
}

// Derive a function signature from a `.ivp` model using diff-grok's parser. Builds the
// `//input:` / `//output:` lines directly from the parsed object (input annotations are stored
// verbatim on it), matching the runnable script's names exactly — no codegen needed.
function preparseIvpModel(parser, text) {
  const {getIVP, MAX_LINE_CHART, STAGE_COL_NAME} = parser;
  let ivp;
  try {
    ivp = getIVP(text);
  } catch {
    return null;
  }
  if (!ivp || !ivp.name) return null;

  // Input names use the script-form names (arg bounds/step `_`-prefixed, loop count `_count`) so the
  // run path, the diff-grok fitting/SA pipeline, and `propagateChoice` lookups all agree on one set
  // of names. `scriptKey` mirrors `name`; it is kept only for the value-forwarding map below.
  // Numeric IVP inputs are always required: mark them `nullable: false` so an emptied field fails
  // form validation instead of running with a null (which throws in the solver).
  const withNullable = (annot) => {
    const a = (annot ?? '').trim();
    if (!a.startsWith('{')) return '{nullable: false}';     // no options block
    if (/[{;\s]nullable\s*:/.test(a)) return a;              // already declares nullable (any value)
    if (/^\{\s*\}$/.test(a)) return '{nullable: false}';     // empty {}
    return `{nullable: false; ${a.slice(1).trimStart()}`;    // everything after '{' untouched
  };

  const mk = (type, name, scriptKey, input) =>
    ({tsType: 'number', name, scriptKey,
      annotation: `//input: ${type} ${name} = ${input.value} ${withNullable(input.annot)}`.trim()});

  const inputs = [];
  const a = ivp.arg.name;
  if (ivp.loop) inputs.push(mk('int', '_count', '_count', ivp.loop.count));
  inputs.push(mk('double', `_${a}0`, `_${a}0`, ivp.arg.initial));
  inputs.push(mk('double', `_${a}1`, `_${a}1`, ivp.arg.final));
  inputs.push(mk('double', '_h', '_h', ivp.arg.step));
  for (const [k, v] of ivp.inits) inputs.push(mk('double', k, k, v));
  if (ivp.params) for (const [k, v] of ivp.params) inputs.push(mk('double', k, k, v));

  const cols = ivp.outputs ? ivp.outputs.size : ivp.inits.size;
  const multiAxis = cols > MAX_LINE_CHART - 1 ? 'true' : 'false';
  const segments = ivp.updates ? ` segmentColumnName: "${STAGE_COL_NAME}",` : '';
  // The 'DiffStudio Facet' viewer (one line chart per variable) is appended last; it resolves
  // only when the DiffStudio package is installed, so Grid / Line chart remain the always-available
  // viewers for packages that ship `.ivp` models without depending on DiffStudio.
  const facetSegment = ivp.updates ? `, segmentColumnName: "${STAGE_COL_NAME}"` : '';
  const viewer = `Grid(block: 100) | Line chart(block: 100, multiAxis: "${multiAxis}",${segments} ` +
    `multiAxisLegendPosition: "RightCenter", autoLayout: "false", showAggrSelectors: "false")` +
    ` | DiffStudio Facet(block: 100${facetSegment})`;
  const outputAnnotation = `//output: dataframe df {caption: ${ivp.name}; viewer: ${viewer}}`;

  const metas = {runOnOpen: 'true', runOnInput: 'true', features: '{"sens-analysis": true, "fitting": true}'};
  let isModel = false;
  for (const meta of ivp.metas || []) {
    const m = /^(?:\/\/)?meta\.([\w-]+):\s*(.*)$/.exec(String(meta).trim());
    if (!m) continue;
    if (m[1] === 'role') {
      if (m[2].split(',').map((s) => s.trim()).includes('model')) isModel = true;
    } else if (m[1] !== 'solver') metas[m[1]] = m[2];
  }
  // Convert a `#meta.inputs` lookup into a real `propagateChoice: all` string input. RFV renders
  // this natively (it ignores `meta.inputs`); selecting a row fills the matching inputs by name.
  // Name / brace parsing mirrors DiffStudio's getLookupsInfo (utils.ts): the name is the text
  // before the first `{` with spaces stripped, and the options block ends at the first `}`.
  if (ivp.inputsLookup) {
    const raw = ivp.inputsLookup.trim();
    const open = raw.indexOf('{');
    const close = raw.indexOf('}');
    if (open >= 0 && close >= 0) {
      const name = raw.slice(0, open).replaceAll(' ', '');
      inputs.unshift({
        tsType: 'string', name, scriptKey: null,
        annotation: `//input: string ${raw.slice(0, close)}; propagateChoice: all${raw.slice(close)}`.trim(),
      });
    }
  }

  return {name: ivp.name, description: ivp.descr || undefined, isModel, inputs, outputAnnotation, metas};
}

const baseImport = 'import * as DG from \'datagrok-api/dg\';\n';

function normEol(s) {
  return s.replace(/\r\n/g, '\n');
}
// eslint-disable-next-line max-len
const annotationForGeneratedFile = `/**\nThis file is auto-generated by the webpack command.\nIf you notice any changes, please push them to the repository.\nDo not edit this file manually.\n*/\n`;

class FuncGeneratorPlugin {
  constructor(options = {outputPath: './src/package.g.ts'}) {
    this.options = options;
  }

  apply(compiler) {
    const srcDirPath = path.join(compiler.context, 'src');
    let packageFilePath = path.join(srcDirPath, 'package.ts');
    if (!fs.existsSync(packageFilePath))
      packageFilePath = path.join(srcDirPath, 'package.js');
    const tsFiles = this._getTsFiles(srcDirPath);
    const genImports = [];
    const genExports = [];

    compiler.hooks.compilation.tap('FuncGeneratorPlugin', (_compilation) => {
      this._clearGeneratedFile();

      for (const file of tsFiles) {
        const content = fs.readFileSync(file, 'utf-8');
        if (!content.includes('@grok.decorators.')) continue;

        if (!content) continue;
        const ast = tsParser.parse(content, {
          sourceType: 'module',
          plugins: [
            ['decorators', {decoratorsBeforeExport: true}],
            'classProperties',
            'typescript',
          ],
        });
        const functions = [];
        const imports = new Set();
        this._walk(ast, (node, parentClass) => {
          const decorators = node.decorators;
          if (!decorators || decorators.length === 0) return;

          if (node?.type === 'ClassDeclaration') {
            this._addNodeData(
              node,
              file,
              srcDirPath,
              functions,
              imports,
              genImports,
              genExports,
            );
          }

          if (node?.type === 'MethodDefinition') {
            this._addNodeData(
              node,
              file,
              srcDirPath,
              functions,
              imports,
              genImports,
              genExports,
              parentClass,
            );
          }
        });
        this._insertImports([...imports]);
        fs.appendFileSync(this.options.outputPath, normEol(functions.join('')), 'utf-8');
      }
      this._generateIvpModels(compiler.context);
      this._checkPackageFileForDecoratorsExport(packageFilePath);
      // Uncommment to add obvious import/export
      // this._writeToPackageFile(packageFilePath, genImports, genExports);
    });
    api({_: ['api']});
  }

  _addNodeData(
    node,
    file,
    srcDirPath,
    functions,
    imports,
    genImports,
    genExports,
    parent = undefined,
  ) {
    if (
      !node.decorators ||
      !node.decorators.length ||
      node.decorators?.length === 0
    )
      return;

    function modifyImportPath(dirPath, filePath) {
      const relativePath = path.relative(dirPath, filePath);
      return `./${relativePath
        .slice(0, relativePath.length - 3)
        .replace(/\\/g, '/')}`;
    }

    const decorator = node.decorators[0];
    const exp = decorator.expression;
    const name = exp.callee?.property?.name || exp.callee?.name;
    const identifierName = node.id?.name || node.key?.name;
    const className = parent?.id?.name || parent?.key?.name;

    if (!name) return;

    const decoratorOptions = this._readDecoratorOptions(
      exp.arguments[0]?.properties ?? [],
    );

    decoratorOptions.set('tags', [
      ...(reservedDecorators[name]['metadata']['tags'] ?? []),
      ...(decoratorOptions.get('tags') ?? []),
    ]);
    
    const role = reservedDecorators[name]['metadata']['role'];
    if (role?.length > 0) {
      const camelRole = toCamelCase(role);

      if (!decoratorOptions.has('meta')) 
        decoratorOptions.set('meta', {role: camelRole});
      else {
        const meta = decoratorOptions.get('meta');
        meta['role'] = meta['role'] ? `${meta['role']},${camelRole}` : camelRole;
      }
    }

    const functionParams =
      node?.type === 'MethodDefinition' ? this._readMethodParamas(node) : [];
    const annotationByReturnObj = node?.type === 'MethodDefinition' ? this._readOutputsFromReturnType(node) : undefined;
    
    const isMethodAsync = this._isMethodAsync(node);
    const importString = generateImport(
      node?.type === 'MethodDefinition' ? className : identifierName,
      modifyImportPath(path.dirname(this.options.outputPath), file),
    );
    const metadataCopy = {...reservedDecorators[name]['metadata']};
    delete metadataCopy.role;
    
    imports.add(importString);
    const funcName = `${
      node?.type === 'MethodDefinition' ? '' : '_'
    }${identifierName}`;
    const funcAnnotaionOptions = {
      ...{name: funcName},
      ...metadataCopy,
      ...(annotationByReturnObj ?
        {outputs: annotationByReturnObj ?? []} :
        {}),
      ...Object.fromEntries(decoratorOptions),
      ...{inputs: functionParams},
      ...{isAsync: isMethodAsync},
    };

    let actualType = undefined;
    if (annotationByReturnObj?.length > 1)
      actualType = 'any';
    else if (annotationByReturnObj?.length === 1 && annotationByReturnObj[0]?.actualType)
      actualType = annotationByReturnObj[0].actualType;
    else if (funcAnnotaionOptions.outputs?.length === 0)
      actualType = 'void';

    // if (!funcAnnotaionOptions.name) funcAnnotaionOptions.name = identifierName;
    function containsAnything(funcAnnotaionOptions) {
      let hasValues = false;
      const arrays = ['tags', 'inputs', 'outputs'];
      for (const option of Object.keys(funcAnnotaionOptions)) {
        if (arrays.includes(option)) {
          if (funcAnnotaionOptions[option].length > 0)
            hasValues = true;
        } else if (option != 'isAsync' && option != 'name')
          hasValues = true;
      }
      return hasValues;
    }
    if (funcAnnotaionOptions.name === funcName && containsAnything(funcAnnotaionOptions))
      delete funcAnnotaionOptions.name;
    functions.push(
      reservedDecorators[name]['genFunc'](
        getFuncAnnotation(funcAnnotaionOptions),
        identifierName,
        '\n',
        className ?? '',
        functionParams,
        actualType,
        funcAnnotaionOptions.isAsync ?? false,
      ),
    );
    genImports.push(
      generateImport(
        funcName,
        modifyImportPath(srcDirPath, this.options.outputPath),
      ),
    );
    genExports.push(generateExport(funcName));
  }

  _isMethodAsync(node) {
    let result = false;
    if (node.type === 'MethodDefinition') result = node.value.async;
    return result;
  }

  _readImports(content) {
    const importRegex = /(import(?:[\s\w{},*]+from\s*)?['']([^'']+)[''];)/g;
    const results = [];

    let match;
    while ((match = importRegex.exec(content)) !== null) 
      results.push(`${match[1]}\n`);
    
    return results;
  }

  _readDecoratorOptions(properties, readForParams= false) {
    const resultMap = new Map();
    let resultObj = undefined;
    // for (const prop of properties)
    //   resultMap.set(prop.key.name ?? prop.key.value, this._evalLiteral(prop.value));
    const optionsToAdd = new Map();
    for (const prop of properties) {
      const key = decoratorOptionToAnnotation.get(prop.key.name ?? prop.key.value ?? '') ?? prop.key.name ?? prop.key.value;
      if (!readForParams && key === 'result')
        resultObj = this._evalLiteral(prop.value);
      else if (readForParams && inputOptionsNames.includes( prop.key.name ?? prop.key.value)) 
        optionsToAdd.set(key, this._evalLiteral(prop.value));
      else
        resultMap.set(key, this._evalLiteral(prop.value));
    }
    
    for (const [key, value] of optionsToAdd)
      resultMap.set(key, value);
    if (resultObj)
      resultMap.set('outputs', [resultObj]);
    
    return resultMap;
  }

  _evalLiteral(node) {
    if (!node) return null;
    switch (node.type) {
    case 'Literal':
      return node.value;

    case 'ArrayExpression':
      return node.elements.map((el) => this._evalLiteral(el));

    case 'MemberExpression':
      return dgAnnotationTypes[node?.property?.name ?? ''] ?? 'dynamic';

    case 'ObjectExpression':
      return Object.fromEntries(
        node.properties.map((p) => {
          return [p.key.name || p.key.value, this._evalLiteral(p.value)];
        }),
      );

    default:
      return '';
    }
  }

  _readMethodParamas(node) {
    const params = node?.value?.params?.map((param) => {
      let baseParam =
        param.type === 'TSNonNullExpression' ? param.expression : param;
      const options =
        param.decorators?.length > 0 ?
          Object.fromEntries(
            this._readDecoratorOptions(
              param.decorators[0]?.expression?.arguments[0].properties, true),
          ) :
          undefined;
        
      // Commented code finds value by default of function's variable
      // let defaultValue = undefined;
      if (baseParam.type === 'AssignmentPattern') {
        // if (baseParam?.right?.type  === 'Literal')
        //   defaultValue = baseParam?.right?.raw;
        baseParam = baseParam?.left;
      }
      const optional = param.optional;
      if (baseParam.type === 'RestElement' || baseParam.type === 'Identifier') {
        let name =
          baseParam.type === 'RestElement' ?
            `...${baseParam.argument.name}` :
            baseParam.name;
        if (baseParam?.argument?.typeAnnotation) {
          name +=
            ': ' +
            generate(baseParam.argument.typeAnnotation.typeAnnotation).code;
        }

        let type = '';
        if (baseParam?.typeAnnotation?.typeAnnotation)
          type = generate(baseParam.typeAnnotation.typeAnnotation).code;
        else type = 'any';

        const params =
          baseParam.typeAnnotation.typeAnnotation.typeArguments?.params;
        if (type !== 'any' && params && params.length > 0)
          type += `<${params.map((e) => e.typeName?.name ?? 'any').join(',')}>`;
        return {name: name, type: type, options: options, optional: optional};
      }
      // Commented code belove sets more strong types for ObjectPatterns and ArrayPatterns
      // else if (baseParam.type === 'ObjectPattern' || baseParam.type === 'ArrayPattern') {
      //   let name = '';
      //   if (baseParam.type === 'ObjectPattern') {
      //     const properties = baseParam.properties.map(prop => {
      //       if (prop.type === 'Property' && prop.key.type === 'Identifier')
      //         return prop.key.name;
      //       else if (prop.type === 'RestElement' && prop.argument.type === 'Identifier')
      //         return `...${prop.argument.name}`;
      //       else
      //         return generate(prop).code;
      //     });
      //     name = `{ ${properties.join(', ')} }`;
      //   } else {
      //     const elements = baseParam.elements.map(elem => {
      //       if (elem) {
      //         if (elem.type === 'Identifier')
      //           return elem.name;
      //         else
      //           return generate(elem).code;
      //       } else return '';
      //     });
      //     name = `[${elements.join(', ')}]`;
      //   }

      //   let type = '';
      //   if (baseParam.typeAnnotation)
      //     type = generate(baseParam.typeAnnotation.typeAnnotation).code;
      //   else
      //     type = 'any';

      //   return { name: name, type: type, options: options };
      // }

      return {name: 'value', type: 'any', options: undefined, optional: optional};
    });
    return params;
  }

  _getTsFiles(root) {
    const tsFiles = [];
    const extPattern = /\.tsx?$/;
    const excludedFiles = ['package-test.ts', 'package.g.ts'];

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
    if (!node || typeof node !== 'object') return;

    visitor(node, parent);

    for (const key in node) {
      const value = node[key];
      if (Array.isArray(value)) {
        value.forEach((child) => {
          if (child && typeof child.type === 'string') {
            this._walk(
              child,
              visitor,
              node.type === 'ClassDeclaration' ? node : parent,
            );
          }
        });
      } else if (value && typeof value.type === 'string') {
        this._walk(
          value,
          visitor,
          node.type === 'ClassDeclaration' ? node : parent,
        );
      }
    }
  }

  _readReturnType(annotation) {
    let resultType = 'void';
    let nodeAnnotation = annotation;
    let isArray = false; 
    let isNullable = false; 
    if (nodeAnnotation?.type === 'TSUnionType' && 
      nodeAnnotation?.types?.length === 2 && 
      nodeAnnotation?.types?.some((e)=> e?.type === 'TSNullKeyword' || e?.type === 'TSVoidKeyword'|| e?.type === 'TSUndefinedKeyword')) {
      nodeAnnotation = 
        nodeAnnotation.types.filter((e)=> e.type !== 'TSNullKeyword' || e?.type === 'TSVoidKeyword' || e?.type === 'TSUndefinedKeyword')[0];
      isNullable = true;
    }
    

    if (
      nodeAnnotation &&
      nodeAnnotation.type !== 'TSUnionType' &&
      nodeAnnotation.type !== 'TSIntersectionType'
    ) {
      if (nodeAnnotation.typeName || nodeAnnotation.type === 'TSTypeReference') {
        resultType =
          nodeAnnotation.typeName?.right?.name ?? nodeAnnotation.typeName?.name;
      } else if (nodeAnnotation.type !== 'TSArrayType') 
        resultType = this._getTypeNameFromNode(nodeAnnotation);
      else if (nodeAnnotation.elementType.type !== 'TSTypeReference') {
        isArray = true;
        resultType = this._getTypeNameFromNode(nodeAnnotation?.elementType);
      } else {
        isArray = true;
        resultType =
          nodeAnnotation?.elementType?.typeName?.name ||
          nodeAnnotation?.elementType?.typeName?.right?.name;
      } 
    }
    let anotationType = typesToAnnotation[resultType];
    resultType = isNullable? 'any' : resultType;
    if (isArray && anotationType) anotationType = `list<${anotationType}>`;
    return [anotationType ?? 'dynamic', `${resultType}${isArray? '[]': ''}` ?? 'any'];
  }

  _readOutputsFromReturnType(node) {
    let results = [];
    let annotation = node.value?.returnType?.typeAnnotation;

    if (node?.type === 'ClassDeclaration')
      return [];
    
    if (annotation?.typeName?.name === 'Promise') {
      const argumnets = annotation.typeArguments?.params;
      if (argumnets && argumnets.length === 1) 
        annotation = argumnets[0];
      else annotation = {};
    }

    if (annotation?.type === 'TSTypeLiteral') 
      results = this._readOutputsFromReturnTypeObject(annotation);
    else {
      const [annotationType, resultType] = this._readReturnType(annotation);
      if (annotationType !== 'void') {
        results.push({
          name: 'result', 
          type: annotationType,
          actualType: resultType,
        });
      }
    }
    return results;
  }

  _readOutputsFromReturnTypeObject(node) {
    let i = 0;
    const results = [];
    for (const member of node.members) {
      const [annotationType, resultType] = this._readReturnType(annotation);
      results.push({
        name: member?.key?.name ?? `result${i}`,
        type: annotationType,
        actualType: resultType,
      });
      i++;
    }
    return results;
  }

  _getTypeNameFromNode(typeNode) {
    if (typeNode.type === 'TSTypeReference') 
      return typeNode.typeName.name;
    else if (typeNode.type === 'TSVoidKeyword') 
      return 'void';
    else if (typeNode.type === 'TSUndefinedKeyword ') 
      return 'undefined';
    else if (typeNode.type === 'TSNumberKeyword') 
      return 'number';
    else if (typeNode.type === 'TSStringKeyword') 
      return 'string';
    else if (typeNode.type === 'TSBooleanKeyword') 
      return 'boolean';
    else 
      return typeNode.type;
    
  }

  // Scan the package's `files/` tree for `.ivp` models. Each one with a `#meta.role: model`
  // line is preparsed into a function signature and emitted as a Rich-Function-View model
  // wrapper that links back to the deployed `.ivp` (multithreaded fitting reads that link).
  // Parse failures are warned, never thrown, so a malformed `.ivp` cannot break the build.
  _generateIvpModels(context) {
    const filesDir = path.join(context, 'files');
    if (!fs.existsSync(filesDir)) return;

    const parser = loadIvpParser();
    if (!parser) return;

    const ivpFiles = fs.readdirSync(filesDir, {recursive: true})
      .filter((f) => typeof f === 'string' && f.endsWith('.ivp'))
      .map((f) => path.join(filesDir, f));

    const namespace = path.basename(context);
    const blocks = [];
    for (const file of ivpFiles) {
      const rel = path.relative(context, file).replace(/\\/g, '/');
      let res;
      try {
        res = preparseIvpModel(parser, fs.readFileSync(file, 'utf-8'));
      } catch (e) {
        console.warn(`[func-gen] skipped ivp model ${rel}: ${e.message}`);
        continue;
      }
      if (!res) {
        console.warn(`[func-gen] skipped ivp model ${rel}: missing #name`);
        continue;
      }
      if (!res.isModel) continue;

      const modelPath = `System:AppData/${namespace}/${path.relative(filesDir, file).replace(/\\/g, '/')}`;
      const fnName = 'ivpModel_' + res.name.replace(/[^A-Za-z0-9]/g, '_');
      const argList = res.inputs.map((i) => `${i.name}: ${i.tsType}`).join(', ');
      const argObj = res.inputs
        .filter((i) => i.scriptKey)
        .map((i) => i.scriptKey === i.name ? i.name : `${i.scriptKey}: ${i.name}`)
        .join(', ');
      // res.metas already merges the shared Diff Studio defaults with the model's own #meta.*.
      const meta = {
        ...res.metas,
        role: 'model',
        diffStudioModel: modelPath,
      };
      // Emit name / description / editor / meta through the canonical annotation generator
      // (handles meta key translation); keep the IVP-derived //input / //output lines verbatim
      // so they match the runnable script exactly.
      const header = getFuncAnnotation({
        name: res.name,
        description: res.description,
        editor: 'Compute2:RichFunctionViewEditor',
        meta,
        inputs: [],
        outputs: [],
      }).trimEnd();
      const annotations = [header, ...res.inputs.map((i) => i.annotation), res.outputAnnotation].join('\n');
      blocks.push(`${annotations}\nexport async function ${fnName}(${argList}): Promise<DG.DataFrame> {\n` +
        `  return await runDiffStudioModel('${modelPath}', {${argObj}});\n}\n`);
    }

    if (blocks.length === 0) return;
    const importLine = `import {runDiffStudioModel} from './ivp-runtime';\n`;
    let content = fs.readFileSync(this.options.outputPath, 'utf-8');
    if (!content.includes(importLine)) content = content.replace(baseImport, baseImport + importLine);
    content += '\n' + blocks.join('\n');
    fs.writeFileSync(this.options.outputPath, normEol(content), 'utf-8');
  }
  _clearGeneratedFile() {
    fs.writeFileSync(this.options.outputPath, normEol(baseImport));
  }

  _checkPackageFileForDecoratorsExport(packagePath) {
    const content = fs.readFileSync(packagePath, 'utf-8');
    const decoratorsExportRegex = /export\s*\*\s*from\s*'\.\/package\.g';/;
    if (!decoratorsExportRegex.test(content)) {
      console.warn(
        // eslint-disable-next-line max-len
        `\nWARNING: Your package doesn't export package.g.ts file to package.ts \n please add 'export * from './package.g';' to the package.ts file.\n`,
      );
    }
  }

  _writeToPackageFile(filePath, imports, exp) {
    if (imports.length !== exp.length) return;
    let content = fs.readFileSync(filePath, 'utf-8');
    for (let i = 0; i < imports.length; i++) {
      const importStatement = imports[i];
      const exportStatement = exp[i];
      if (!content.includes(importStatement.trim()))
        content = annotationForGeneratedFile + importStatement + content + exportStatement;
    }
    fs.writeFileSync(filePath, content, 'utf-8');
  }

  _insertImports(importArray) {
    if (fs.existsSync(this.options.outputPath)) {
      const content = fs.readFileSync(this.options.outputPath, 'utf-8');
      if (content) importArray.push(content);
      const output = importArray.join('');
      fs.writeFileSync(this.options.outputPath, normEol(output), 'utf-8');
    } else {
      fs.writeFileSync(
        this.options.outputPath,
        normEol(`${annotationForGeneratedFile}${baseImport}\n${importArray.join('')}`),
        'utf-8',
      );
    }
  }
}

module.exports = FuncGeneratorPlugin;
