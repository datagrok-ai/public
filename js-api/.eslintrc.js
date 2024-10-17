// eslint-disable-next-line no-undef
module.exports = {
  'env': {
    'browser': true,
    'es2021': true
  },
  "plugins": [
    "@datagrok-misc/config"
  ],
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
  }
};
