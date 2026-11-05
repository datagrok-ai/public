/* The EMS gallery over u2's static server: served from `libraries/`, so the pages reach u2's built
   `src/` and sheets through the import map in `gallery/index.html`. */
import {resolve} from 'node:path';
import {fileURLToPath} from 'node:url';
import {startGalleryServer as serve} from '../../u2/e2e/gallery-server.mjs';

const LIBRARIES = resolve(fileURLToPath(new URL('../..', import.meta.url)));
const SERVED = ['ems/gallery', 'ems/src', 'ems/css', 'u2/src', 'u2/css', 'u2/node_modules/@floating-ui',
  'u2/node_modules/datagrok-api/dist/src/u2core', 'u2/node_modules/datagrok-api/dist/src/const.js'];

export const startGalleryServer = () => serve(LIBRARIES, SERVED);
