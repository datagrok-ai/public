import fs from 'fs';
import path from 'path';
import walk from 'ignore-walk';
import * as utils from '../utils/utils';
import * as color from '../utils/color-utils';
import { ValidationResult } from '../validators/interfaces';


export function check(args: { [x: string]: string | string[]; }) {
  const nOptions = Object.keys(args).length - 1;
  if (args['_'].length !== 1 || nOptions > 0)
    return false;

  const curDir = process.cwd();
  if (!utils.isPackageDir(curDir)) {
    color.error('File `package.json` not found. Run the command from the package directory');
    return false;
  }

  const files = walk.sync({ ignoreFiles: ['.npmignore', '.gitignore'] })
    .filter((f) => f.endsWith('.js') || f.endsWith('.ts'));

  const webpackConfigPath = path.join(curDir, 'webpack.config.js');
  const isWebpack = fs.existsSync(webpackConfigPath);
  if (isWebpack) {
    const content = fs.readFileSync(webpackConfigPath, { encoding: 'utf-8' });
    const externals = extractExternals(content);
    if (externals)
      checkImportStatements(files, externals).forEach((warning) => color.warn(warning));
  }

  checkFuncSignatures(files).forEach((warning) => color.warn(warning));

  return true;
}

function extractExternals(config: string): {}|null {
  const externalsRegex = /(?<=externals\s*:\s*)(\{[\S\s]*?\})/;
  const match = config.match(externalsRegex);
  if (match) {
    // Replace single quotes and a trailing comma to make a string JSON-like
    const externalStr = match[1].replace(/'/g, '"').replace(/(?<=[\S\s]),(?=\s*\})/, '');
    try {
      const externals = JSON.parse(externalStr);
      return externals;
    } catch(e) {
      return null;
    }
  }
  return null;
}

function checkImportStatements(files: string[], externals: {}): string[] {
  const modules = [];
  for (const key in externals) {
    modules.push(key);
    if (key.includes('/'))
      modules.push(key.split('/', 1)[0]);
  }
  const importRegex = new RegExp(`import\\s+.*(${modules.join('|')}).*(?=\\s+?)`, 'g');
  const validImportRegex = new RegExp(`import\\s+.*(${Object.keys(externals).join('|')})['"]{1}`);
  const warnings: string[] = [];

  function validateImport(file: string, s: string): ValidationResult {
    let value = validImportRegex.test(s);
    let message = value ? '' : 'Pay attention to file ' + file + ': import statement `' +
      s + '` differs from the path given in the webpack config as an external module. ' +
      'It can increase the bundle size.';
    return { value, message };
  }

  for (const file of files) {
    const content = fs.readFileSync(path.join(process.cwd(), file), { encoding: 'utf-8' });
    const matchedImports = content.match(importRegex);
    if (matchedImports) {
      for (const match of matchedImports) {
        const vr = validateImport(file, match);
        if (!vr.value)
          warnings.push(vr.message);
      }
    }
  }

  return warnings;
}

function checkFuncSignatures(files: string[]): string[] {
  const warnings: string[] = [];
  const checkFunctions: { [role: string]: FuncValidator } = {
    semTypeDetector: ({inputs, outputs}: {inputs: FuncParam[], outputs: FuncParam[]}) => {
      let value = true;
      let message = '';

      if (inputs.length !== 1 || inputs[0].type !== 'column') {
        value = false;
        message += 'Semantic type detectors must have one input of type "column"';
      }

      if (outputs.length !== 1 || outputs[0].type !== 'string') {
        value = false;
        message += 'Semantic type detectors must have one output of type "string"';
      }

      return { value, message };
    },
    cellRenderer: ({inputs, outputs}: {inputs: FuncParam[], outputs: FuncParam[]}) => {
      let value = true;
      let message = '';

      if (inputs.length !== 0) {
        value = false;
        message += 'Cell renderer functions should take no arguments';
      }

      if (outputs.length !== 1 || outputs[0].type !== 'grid_cell_renderer') {
        value = false;
        message += 'Cell renderer functions must have one output of type "grid_cell_renderer"';
      }

      return { value, message };
    },
    fileViewer: ({inputs, outputs, tags}: {inputs: FuncParam[], outputs: FuncParam[], tags?: string[]}) => {
      let value = true;
      let message = '';

      if (tags == null || tags.filter((t) => t.startsWith('fileViewer')).length < 2) {
        value = false;
        message += 'File viewers must have at least two special tags: "fileViewer" and "fileViewer-<extension>"';
      }

      if (inputs.length !== 1 || inputs[0].type !== 'file') {
        value = false;
        message += 'File viewers must have one input of type "file"';
      }

      if (outputs.length !== 1 || outputs[0].type !== 'view') {
        value = false;
        message += 'File viewers must have one output of type "view"';
      }

      return { value, message };
    },
    fileExporter: ({description}: {description?: string}) => {
      let value = true;
      let message = '';

      if (description == null || description === '') {
        value = false;
        message += 'File exporters should have a description parameter';
      }

      return { value, message };
    },
    packageSettingsEditor: ({outputs}: {outputs: FuncParam[]}) => {
      let value = true;
      let message = '';

      if (outputs.length === 1 && outputs[0].type === 'widget') {
        value = false;
        message += 'Package settings editors must have one output of type "widget"';
      }

      return { value, message };
    },
  };
  const functionRoles = Object.keys(checkFunctions);

  for (const file of files) {
    const content = fs.readFileSync(path.join(process.cwd(), file), { encoding: 'utf-8' });
    const functions = getFuncMetadata(content);
    for (const f of functions) {
      const roles = functionRoles.filter((role) => f.tags?.includes(role));
      if (roles.length > 1) {
        warnings.push(`File ${file}, function ${f.name}: several function roles are used (${roles.join(', ')})`);
      } else if (roles.length === 1) {
        const vr = checkFunctions[roles[0]](f);
        if (!vr.value)
          warnings.push(`File ${file}, function ${f.name}: ${vr.message}`);
      }
    }
  }

  return warnings;
}

function getFuncMetadata(script: string): FuncMetadata[] {
  //TODO: retrieve function metadata from file
  return [];
}

type FuncParam = {type: string};

type FuncMetadata = {
  name: string,
  inputs: FuncParam[],
  outputs: FuncParam[],
  tags?: string[],
  description?: string,
};

type FuncValidator = ({}: {
  name: string,
  inputs: FuncParam[],
  outputs: FuncParam[],
  tags?: string[],
  description?: string,
}) => ValidationResult;
