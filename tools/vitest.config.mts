import {defineConfig} from 'vitest/config';

export default defineConfig({
  // Babel writes each .js beside its .ts in bin/: tests must load the sources, and one copy of each module.
  resolve: {extensions: ['.ts', '.mts', '.mjs', '.js', '.jsx', '.tsx', '.json']},
  test: {
    projects: [
      {
        extends: true,
        test: {
          name: 'unit',
          environment: 'node',
          include: ['bin/**/*.test.ts'],
          exclude: ['bin/**/*.integration.test.ts'],
        },
      },
      {
        extends: true,
        test: {
          name: 'integration',
          environment: 'node',
          include: ['bin/**/*.integration.test.ts'],
          testTimeout: 30_000,
          hookTimeout: 30_000,
        },
      },
    ],
  },
});
