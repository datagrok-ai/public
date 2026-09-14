// eslint-disable-next-line no-undef
module.exports = {
  'ignorePatterns': ['**/*.d.ts', 'src/api/*.g.ts', 'src/interfaces/*.ts', 'src/datagrok/**', 'node_modules/**'],
  'env': {
    'browser': true,
    'es2021': true
  },
  'parserOptions': {
    'ecmaVersion': 12,
    'sourceType': 'module'
  },
  'rules': {
    'indent': [
      'error',
      2
    ],
    'linebreak-style': [
      'error',
      'windows'
    ],
    'quotes': [
      'error',
      'single'
    ],
    'semi': [
      'error',
      'always'
    ]
  },
  // The TypeScript sources: only the rules that guard the API conventions (CONVENTIONS.md), so the
  // gate stays green and meaningful; formatting is left to the editor.
  'overrides': [{
    'files': ['*.ts', 'src/**/*.ts'],
    'excludedFiles': ['*.d.ts', '**/*.d.ts', 'src/api/*.g.ts', 'src/interfaces/*.ts', 'src/datagrok/**'],
    'parser': '@typescript-eslint/parser',
    'plugins': ['@typescript-eslint'],
    'rules': {
      'indent': 'off',
      'linebreak-style': 'off',
      'quotes': 'off',
      'semi': 'off',
      'no-debugger': 'error',
      '@typescript-eslint/ban-types': ['error', {
        'extendDefaults': false,
        'types': {
          'Function': {'message': 'Declare the real signature, e.g. `() => void` or `(e: MouseEvent) => void`.'},
        },
      }],
    },
  }],
};
