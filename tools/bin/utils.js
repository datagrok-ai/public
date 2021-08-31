const fs = require('fs');

exports.isEmpty = (dir) => fs.readdirSync(dir).length === 0;

exports.kebabToCamelCase = (s) => {
  s = s.replace(/-./g, x => x.toUpperCase()[1]);
  return s[0].toUpperCase() + s.slice(1);
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
};

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
