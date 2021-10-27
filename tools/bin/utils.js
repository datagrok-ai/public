const fs = require('fs');

exports.isEmpty = (dir) => fs.readdirSync(dir).length === 0;

exports.kebabToCamelCase = (s, firstUpper = true) => {
  s = s.replace(/-./g, x => x.toUpperCase()[1]);
  return (firstUpper ? s[0].toUpperCase() : s[0].toLowerCase()) + s.slice(1);
}

exports.spaceToCamelCase = (s, firstUpper = true) => {
  s = s.replace(/\s+./g, x => x[x.length - 1].toUpperCase());
  return (firstUpper ? s[0].toUpperCase() : s[0].toLowerCase()) + s.slice(1);
}

exports.wordsToCamelCase = (s, firstUpper = true) => {
  const m = s.match(/^[A-Z]+/); // only abbreviation handling
  const lastIndex = m ? m[0].length : 1;
  return s.slice(0, lastIndex).toLowerCase() + s.slice(lastIndex);
}

exports.camelCaseToKebab = (s) => {
  return s.replace(/[A-Z]/g, (char, index) => index == 0 ? char.toLowerCase() : '-'+ char.toLowerCase());
}

exports.mapURL = (conf) => {
  let urls = {};
  for (let server in conf.servers) {
    urls[conf['servers'][server]['url']] = server;
  }
  return urls;
}

exports.replacers = {
  NAME: (s, name) => s.replace(/#{NAME}/g, name),
  NAME_TITLECASE: (s, name) => s.replace(/#{NAME_TITLECASE}/g, name[0].toUpperCase() + name.slice(1).toLowerCase()),
  NAME_LOWERCASE: (s, name) => s.replace(/#{NAME_LOWERCASE}/g, name.toLowerCase()),
  NAME_PREFIX: (s, name) => s.replace(/#{NAME_PREFIX}/g, name.slice(0, 3)),
  PACKAGE_DETECTORS_NAME: (s, name) => s.replace(/#{PACKAGE_DETECTORS_NAME}/g, this.kebabToCamelCase(name)),
  PACKAGE_NAMESPACE: (s, name) => s.replace(/#{PACKAGE_NAMESPACE}/g, this.kebabToCamelCase(name)),
  FUNC_NAME: (s, name) => s.replace(/#{FUNC_NAME}/g, name.includes('-') ? this.kebabToCamelCase(name)
    : name.includes(' ') ? this.spaceToCamelCase(name) : name[0].toUpperCase() + name.slice(1)),
  FUNC_NAME_LOWERCASE: (s, name) => s.replace(/#{FUNC_NAME_LOWERCASE}/g, name.includes('-') ?
    this.kebabToCamelCase(name, false) : name.includes(' ') ?
    this.spaceToCamelCase(name, false) : this.wordsToCamelCase(name, false)),
  PARAMS_OBJECT: (s, params) => s.replace(/#{PARAMS_OBJECT}/g, params.length ? `{ ${params.map((p) => p.name).join(', ')} }` : `{}`),
  OUTPUT_TYPE: (s, type) => s.replace(/#{OUTPUT_TYPE}/g, type),
  TYPED_PARAMS: (s, params) => s.replace(/#{TYPED_PARAMS}/g, params.map((p) => `${p.name}: ${p.type}`).join(', '))
};

exports.TemplateBuilder = class TemplateBuilder {
  constructor(template) {
    this.template = template;
  }

  replace(pattern, value) {
    this.template = exports.replacers[pattern](this.template, value);
    return this;
  }

  build() {
    return this.template;
  }
}

exports.scriptLangExtMap = {
  javascript: 'js',
  julia: 'jl',
  node: 'js',
  octave: 'm',
  python: 'py',
  r: 'R',
};

exports.scriptExtensions = ['.jl', '.m', '.py', '.R'];
exports.checkScriptLocation = (filepath) => {
  if (!filepath.startsWith('scripts/') &&
  this.scriptExtensions.some((ext) => filepath.endsWith(ext))) {
    return false;
  }
  return true;
};

exports.getScriptName = (script, comment = '#') => {
  const regex = new RegExp(`${comment}\\s*name:\\s*(.*)`);
  const match = script.match(regex);
  return match ? match[1]?.trim() : null;
};

exports.dgToTsTypeMap = {
  int: 'number',
  double: 'number',
  bigint: 'bigint',
  bool: 'boolean',
  string: 'string',
  dataframe: 'DG.DataFrame',
  column: 'DG.Column',
  column_list: 'string[]',
  file: 'DG.FileInfo',
};

exports.getScriptOutputType = (script, comment = '#') => {
  const regex = new RegExp(`${comment}\\s*output:\\s?([a-z_]+)\\s*`);
  const match = script.match(regex);
  if (!match) return 'void';
  return this.dgToTsTypeMap[match[1]] || 'any';
};

exports.getScriptInputs = (script, comment = '#') => {
  const regex = new RegExp(`${comment}\\s*input:\\s?([a-z_]+)\\s+(\\w+)`, 'g');
  const inputs = [];
  for (let match of script.matchAll(regex)) {
    const type = this.dgToTsTypeMap[match[1]] || 'any';
    const name = match[2];
    inputs.push({ type, name });
  }
  return inputs;
};

exports.dgImports = `import * as grok from 'datagrok-api/grok';\nimport * as DG from 'datagrok-api/dg';\n\n`;

exports.scriptWrapperTemplate = `export async function #{FUNC_NAME_LOWERCASE}(#{TYPED_PARAMS}): Promise<#{OUTPUT_TYPE}> {
  return await grok.functions.call('#{PACKAGE_NAMESPACE}:#{FUNC_NAME}', #{PARAMS_OBJECT});
}`;
