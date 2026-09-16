import {defineConfig} from 'vitest/config';

export default defineConfig({
  test: {
    projects: [
      {
        test: {
          name: 'unit',
          environment: 'node',
          include: ['bin/**/*.test.ts'],
          exclude: ['bin/**/*.integration.test.ts', 'bin/__tests__/fixtures/**'],
          // a `grok kg build` over a fixture monorepo takes seconds and the box is often loaded; no test sets its own limit
          testTimeout: 120_000,
          hookTimeout: 120_000,
        },
      },
      {
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
