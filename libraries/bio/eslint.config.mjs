import typescriptEslint from '@typescript-eslint/eslint-plugin';
import datagrokConfig from '@datagrok-misc/eslint-plugin-config';
import deprecation from 'eslint-plugin-deprecation';
import globals from 'globals';
import tsParser from '@typescript-eslint/parser';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import js from '@eslint/js';
import {FlatCompat} from '@eslint/eslintrc';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);
const compat = new FlatCompat({
  baseDirectory: __dirname,
  recommendedConfig: js.configs.recommended,
  allConfig: js.configs.all,
});

// console.debug(JSON.stringify(compat.extends('google'), undefined, 2));

const globalsBrowserFixed = Object.assign({}, globals.browser, {
  AudioWorkletGlobalScope: globals.browser['AudioWorkletGlobalScope '],
});
delete globalsBrowserFixed['AudioWorkletGlobalScope '];

export default [{
  ignores: ['src/**/*.d.ts*'],
}, ...compat.extends('google'), {
  plugins: {
    '@typescript-eslint': typescriptEslint,
    '@datagrok-misc/config': datagrokConfig,
    deprecation,
  },

  languageOptions: {
    globals: {...globalsBrowserFixed},
    parser: tsParser,
    parserOptions: {project: 'tsconfig.json'},
    ecmaVersion: 12,
    sourceType: 'module',
  },

  rules: {
    'deprecation/deprecation': 'error',
    'indent': ['error', 2, {SwitchCase: 0}],
    'max-len': ['error', 120],
    'no-throw-literal': 'error',
    'no-unused-vars': 'off',
    '@typescript-eslint/no-unused-vars': ['warn', {varsIgnorePattern: '^(_|ui$|grok$|DG$)', argsIgnorePattern: '^_'}],
    'require-jsdoc': 'off',
    'valid-jsdoc': 'warn',
    'spaced-comment': 'off',
    'linebreak-style': 'off',
    'curly': ['error', 'multi-or-nest'],
    'brace-style': ['error', '1tbs', {allowSingleLine: true}],
    'block-spacing': [2, 'always'],
    'comma-dangle': 'off',
    'comma-spacing': 'error',
    'comma-style': 'error',
    'guard-for-in': 'off',
    'space-infix-ops': 'error',
  },

  files: ['**/*.ts'],
}];
