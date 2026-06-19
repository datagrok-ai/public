import js from '@eslint/js';
import tseslint from 'typescript-eslint';
import stylistic from '@stylistic/eslint-plugin';
import globals from 'globals';

// Flat config (ESLint 9). Replaces the legacy .eslintrc.json + eslint-config-google.
// Stylistic rules that used to live in ESLint core now come from @stylistic.
export default tseslint.config(
  {
    ignores: [
      'bin/**/*.js',
      'bin/**/*.js.map',
      'node_modules/**',
      'package-template/**',
      'entity-template/**',
      'script-template/**',
    ],
  },
  js.configs.recommended,
  ...tseslint.configs.recommended,
  {
    files: ['bin/**/*.ts'],
    plugins: {'@stylistic': stylistic},
    languageOptions: {
      ecmaVersion: 2022,
      sourceType: 'module',
      globals: {...globals.browser, ...globals.node},
    },
    rules: {
      '@stylistic/no-trailing-spaces': 'off',
      '@stylistic/indent': ['error', 2],
      '@stylistic/max-len': ['error', 140],
      '@stylistic/padded-blocks': 'off',
      '@stylistic/spaced-comment': 'off',
      '@stylistic/linebreak-style': 'off',
      'guard-for-in': 'off',
      'curly': ['error', 'multi-or-nest'],
      '@stylistic/brace-style': ['error', '1tbs', {allowSingleLine: true}],
      '@stylistic/block-spacing': ['error', 'always'],
    },
  },
);
