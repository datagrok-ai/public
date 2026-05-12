module.exports = {
  preset: 'ts-jest',
  testEnvironment: 'node',
  roots: ['<rootDir>/__tests__'],
  testMatch: ['**/__tests__/**/*.test.ts'],
  testPathIgnorePatterns: ['/node_modules/', '/dist/', '/fixtures/'],
  transform: {
    '^.+\\.ts$': ['ts-jest', {tsconfig: {target: 'es2020', module: 'commonjs', strict: true, esModuleInterop: true}}],
  },
};
