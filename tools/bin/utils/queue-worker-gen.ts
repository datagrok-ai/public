import fs from 'fs';
import path from 'path';
import * as color from './color-utils';

// Matches `//meta.queue: true` in package.ts — functions executed as celery tasks
// in the Node worker (see help/develop/how-to/packages/js-server-functions.md).
const queueTrueRegex = /^\s*\/\/\s*meta\.queue\s*:\s*true\s*$/m;

const QUEUE_DIR = 'queue';

const TEMPLATE_NODE_DOCKER = `FROM datagrok/celery_node_worker:bleeding-edge
EXPOSE 8000
`;

/** Generates dockerfiles/queue/ for packages with \`meta.queue: true\` functions so the
 *  standard docker build/push/image.json flow covers the auto-created Node worker
 *  container (mirror of generateCeleryArtifacts for python/). The worker fetches the
 *  package bundle at runtime, so the image is just the stock worker base. */
export function generateQueueArtifacts(packageDir: string): boolean {
  const candidates = ['package.ts', 'package.g.ts']
    .map((f) => path.join(packageDir, 'src', f))
    .filter((f) => fs.existsSync(f));
  const hasQueueFuncs = candidates.some((f) => queueTrueRegex.test(fs.readFileSync(f, 'utf-8')));
  if (!hasQueueFuncs)
    return false;

  const dockerfilesDir = path.join(packageDir, 'dockerfiles', QUEUE_DIR);
  const dockerfilePath = path.join(dockerfilesDir, 'Dockerfile');
  fs.mkdirSync(dockerfilesDir, {recursive: true});
  if (!fs.existsSync(dockerfilePath))
    fs.writeFileSync(dockerfilePath, TEMPLATE_NODE_DOCKER);

  // Resource config lives in the reserved queue/ package folder
  const containerJsonSrc = path.join(packageDir, QUEUE_DIR, 'container.json');
  const containerJsonDest = path.join(dockerfilesDir, 'container.json');
  if (fs.existsSync(containerJsonSrc) && !fs.existsSync(containerJsonDest))
    fs.copyFileSync(containerJsonSrc, containerJsonDest);

  color.log(`Generated Node worker Docker artifacts in dockerfiles/${QUEUE_DIR}/`);
  return true;
}
