/// The repositories a build reads besides the monorepo and its public submodule: the marketing site (github
/// datagrok-ai/landing), an external root whose paths carry the `landing:` prefix (schema.yaml `Path`). One
/// place turns a prefixed path into a file on disk, a delivery url and a git root, for check, the extractors,
/// the hashing and the browser.
import * as fs from 'fs';
import * as path from 'path';

export const LANDING_PREFIX = 'landing:';
const LANDING_FALLBACK = 'C:/dg/landing';
export const SITE = 'https://datagrok.ai';
const HELP_DIR = 'public/help/';
const DOCUSAURUS_STATIC = 'public/docusaurus/static/';
const WEB = `${LANDING_PREFIX}web/`;
/** The pages nginx serves by their file name; every other `.html` is reached without the extension. */
const KEPT_HTML = ['login.html', 'help.html', 'invite.html', 'js-api.html'];

export interface Roots {
  repoRoot: string;
  landingDir?: string;
}

/** `--landing <dir>`, else `<repo>/../landing`, else the dev-box clone; `--landing false` (the tests) reads no site at all.
 * A folder is the site only when it has `web/`. */
export function landingRoot(repoRoot: string, flag?: unknown): string | undefined {
  if (flag === false || flag === 'false' || flag === 'none') return undefined;
  if (flag !== undefined) return path.resolve(String(flag));
  return [path.resolve(repoRoot, '..', 'landing'), LANDING_FALLBACK].find((d) => fs.existsSync(path.join(d, 'web')));
}

/** Where a repo path (posix, `landing:` prefixed for the site) is on disk; undefined for a site path when no site was given. */
export function localPath(roots: Roots, file: string): string | undefined {
  if (!file.startsWith(LANDING_PREFIX)) return path.join(roots.repoRoot, file);
  return roots.landingDir === undefined ? undefined : path.join(roots.landingDir, file.slice(LANDING_PREFIX.length));
}

export function existsAt(roots: Roots, file: string): boolean {
  const p = localPath(roots, file);
  return p !== undefined && fs.existsSync(p) && fs.statSync(p).isFile();
}

/** The public tree: `public/` and the site. */
export function isPublicPath(p: string): boolean {
  return p.startsWith('public/') || p.startsWith(LANDING_PREFIX);
}

/** Where a reader sees a committed file: the docs site serves help and its static folder under /help, the site its web/ folder at the root. */
export function deliveryUrl(file: string): string | undefined {
  if (file.startsWith(HELP_DIR)) return `${SITE}/help/${file.slice(HELP_DIR.length)}`;
  if (file.startsWith(DOCUSAURUS_STATIC)) return `${SITE}/help/${file.slice(DOCUSAURUS_STATIC.length)}`;
  if (file.startsWith(WEB)) return `${SITE}/${file.slice(WEB.length)}`;
  return undefined;
}

/** The address the site's routing gives a page (`local/nginx.conf`): `.html` dropped, `index.html` the root, a few pages by name. */
export function pageUrl(file: string): string | undefined {
  if (!file.startsWith(WEB)) return undefined;
  const rel = file.slice(WEB.length);
  if (rel === 'index.html') return `${SITE}/`;
  if (KEPT_HTML.includes(rel) || /^[45]0[0234]\.html$/.test(rel)) return `${SITE}/${rel}`;
  return `${SITE}/${rel.replace(/\.html$/i, '')}`;
}

/** The git roots of a build, each with the prefix its paths carry in the graph. */
export function gitRoots(roots: Roots): {dir: string, prefix: string}[] {
  const out = [{dir: roots.repoRoot, prefix: ''}, {dir: path.join(roots.repoRoot, 'public'), prefix: 'public/'}];
  if (roots.landingDir !== undefined) out.push({dir: roots.landingDir, prefix: LANDING_PREFIX});
  return out;
}
