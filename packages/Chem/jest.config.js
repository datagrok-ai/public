module.exports = {
  "roots": [
    "<rootDir>/src"
  ],
  "testMatch": [
    "**/__tests__/**/*.+(ts|tsx)",
    "**/?(*.)+(spec|test).+(ts|tsx)"
  ],
  moduleFileExtensions: [
    'ts',
    'js',
  ],
  "transform": {
    "^.+\\.(ts|tsx)$": "ts-jest"
  },
  transformIgnorePatterns: ['^.+\\.js$'],
  testTimeout: 30000
}