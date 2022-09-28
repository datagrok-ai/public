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

  const webpackConfigPath = path.join(curDir, 'webpack.config.js');
  const isWebpack = fs.existsSync(webpackConfigPath);
  if (isWebpack) {
    const content = fs.readFileSync(webpackConfigPath, { encoding: 'utf-8' });
    const externals = extractExternals(content);
    if (externals) {
      const warnings = checkImportStatements(externals);
      if (warnings.length)
        warnings.forEach((warning) => color.warn(warning));
    }
  }

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

function checkImportStatements(externals: {}): string[] {
  const warnings: string[] = [];
  const files = walk.sync({ ignoreFiles: ['.npmignore', '.gitignore'] });

  function validateImport(s: string): ValidationResult {
    return { value: true, message: ''};
  }

  for (const file of files) {
    if (!file.endsWith('.js') && !file.endsWith('.ts'))
      continue;
    const content = fs.readFileSync(path.join(process.cwd(), file), { encoding: 'utf-8' });
    const importRegex = new RegExp(`import\\s+.*(${Object.keys(externals).join('|')})`, 'g');
    const matchedImports = content.match(importRegex);
    if (matchedImports) {
      for (const match of matchedImports) {
        const vr = validateImport(match);
        if (!vr.value)
          warnings.push(vr.message);
      }
    }
  }

  return warnings;
}
