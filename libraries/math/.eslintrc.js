/* eslint-disable max-len */

module.exports = {
  env: {
    browser: true,
    es6: true,
  },
  extends: ['google'],
  parser: '@typescript-eslint/parser',
  parserOptions: {
    ecmaVersion: 12,
    sourceType: 'module',
  },
  plugins: ['@typescript-eslint', 'gpu-rules', 
    "@datagrok-misc/config"
  ],
  rules: {
    'gpu-rules/assure-buffer-destroy': 'error',
    'indent': ['error', 2, {SwitchCase: 1}],
    'max-len': ['error', 120],
    'no-unused-vars': 'off',
    '@typescript-eslint/no-unused-vars': [
      'warn',
      {varsIgnorePattern: '^(_|ui$|grok$|DG$)', argsIgnorePattern: '^_'},
    ],
    'require-jsdoc': 'off',
    'valid-jsdoc': 'off',
    'spaced-comment': 'off',
    'linebreak-style': 'off',
    'curly': ['error', 'multi-or-nest', 'consistent'],
    'brace-style': ['error', '1tbs', {allowSingleLine: true}],
    'block-spacing': [2, 'always'],
    'comma-dangle': [
      'error',
      {
        arrays: 'only-multiline',
        functions: 'never',
        objects: 'only-multiline',
        imports: 'only-multiline',
      },
    ],
    'guard-for-in': 'off',
    'no-restricted-syntax': [
      'error',
      {
        selector:
          'MemberExpression[object.name=\'device\'][property.name=\'destroy\']',
        message:
          'Use of device.destroy is disallowed. device is used in singleton manner. make sure to destroy all buffers on that GPU before returning.',
      },
      {
        selector:
          'MemberExpression[object.property.name=\'gpu\'][property.name=\'requestAdapter\']',
        message:
          'Use of gpu.requestAdapter is disallowed. device is used in singleton manner. use getGPUDevice method instead.',
      },
    ],
  },
};
