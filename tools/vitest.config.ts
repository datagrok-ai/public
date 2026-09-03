import {defineConfig} from 'vitest/config';

// After `npm run build` every bin/**/*.ts has a compiled .js twin next to it, and Vite's
// default extension order would resolve `../commands/x` to that (stale, CommonJS) twin —
// bypassing vi.mock and testing the previous build. Tests always run from the sources.
const resolve = {extensions: ['.ts', '.mts', '.js', '.mjs', '.json']};

export default defineConfig({
  test: {
    projects: [
      {
        resolve,
        test: {
          name: 'unit',
          environment: 'node',
          include: ['bin/**/*.test.ts'],
          exclude: ['bin/**/*.integration.test.ts'],
        },
      },
      {
        resolve,
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
