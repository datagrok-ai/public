import {defineConfig} from '@playwright/test';
import {baseConfig} from '@datagrok-libraries/test/src/playwright/base-config';

// BiostructureViewer package-owned E2E suite. General config lives in the shared base; only run-dir
// specifics belong here.
export default defineConfig({
  ...baseConfig,
  testDir: '.',
  use: {
    ...baseConfig.use,
    launchOptions: {
      ...baseConfig.use?.launchOptions,
      // Mol*/molstar (3D Structure pane) needs a WebGL2 context; the GPU-less CI container has none,
      // so the viewer never mounts (.bsv-container-info-panel absent) though it renders locally. Force
      // software WebGL via ANGLE + SwiftShader so the 3D structure viewer mounts on CI too.
      args: [
        ...(baseConfig.use?.launchOptions?.args ?? []),
        '--use-gl=angle',
        '--use-angle=swiftshader',
        '--enable-unsafe-swiftshader',
        '--ignore-gpu-blocklist',
      ],
    },
  },
});
