module.exports = {
  'roots': [
    '<rootDir>/src',
  ],
  'testMatch': [
    '**/__jest__/**/*.test.+(ts|tsx)',
  ],
  moduleFileExtensions: [
    'ts',
    'js',
  ],
  'transform': {
    '^.+\\.(ts|tsx)$': 'ts-jest',
  },
  transformIgnorePatterns: ['^.+\\.js$'],
  globals: {
    'ts-jest': {
      'tsconfig': {
        'target': 'es6',
        'module': 'es2020',
      },
    },
  },
  reporters: [
    'default',
    [
      './node_modules/jest-html-reporter',
      {
        'includeConsoleLog': true,
      },
    ],
  ],
};
