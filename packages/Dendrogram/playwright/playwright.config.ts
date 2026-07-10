import {defineConfig} from '@playwright/test';
import {baseConfig} from '@datagrok-libraries/test/src/playwright/base-config';

// Dendrogram package-owned E2E suite. General config lives in the shared base; only run-dir
// specifics belong here.
export default defineConfig({
  ...baseConfig,
  testDir: '.',
  use: {
    ...baseConfig.use,
    launchOptions: {
      ...baseConfig.use?.launchOptions,
      // The GPU-less CI container has no hardware WebGL; PhylocanvasGL (regl/WebGL) fails there
      // with "Bad state: No element" though it renders fine locally. Force software WebGL via ANGLE
      // + SwiftShader so the phylogenetic-tree canvas mounts on CI too.
      args: [
        ...(baseConfig.use?.launchOptions?.args ?? []),
        '--use-gl=angle',
        '--use-angle=swiftshader',
        '--enable-unsafe-swiftshader',
      ],
    },
  },
});
