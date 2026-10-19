/* u2's stylesheet, installed once per page. Every package that bundles u2 carries a copy; the
   first to load inserts one <style data-u2-styles>, and a later copy replaces its text only when its
   revision is newer — so the page always runs the newest u2 CSS, whatever order packages load in,
   and a CSS fix lands on a stand as soon as one rebuilt package does. Nothing to import: loading
   any u2 control installs it. A standalone host that links css/*.css itself opts out with
   <meta name="u2-styles" content="linked">. */
import {U2_STYLES} from './styles.g.js';

export interface InstalledStyles {
  revision: number;
  hash: string;
  element: HTMLStyleElement;
}

const SELECTOR = 'style[data-u2-styles]';
let warned = false;
let legacyWarned = false;

export class Styles {
  /** This copy's revision and content hash. */
  static readonly revision: number = U2_STYLES.revision;
  static readonly hash: string = U2_STYLES.hash;

  /** Installs this copy's CSS unless the page already runs the same or a newer revision. */
  static install(css: string = U2_STYLES.css, revision: number = U2_STYLES.revision,
    hash: string = U2_STYLES.hash): InstalledStyles | null {
    const doc = globalThis.document;
    if (doc?.head === undefined || doc.head === null ||
      doc.querySelector('meta[name="u2-styles"][content="linked"]') !== null)
      return null;
    if (!legacyWarned && Styles.legacyCopies(doc) > 0) {
      legacyWarned = true;
      console.warn('u2: a package built before stylesheet revisions injects its own copy of the u2 CSS, which overrides ' +
        `revision ${revision}; republish the packages that bundle u2.`);
    }
    let el = doc.head.querySelector<HTMLStyleElement>(SELECTOR);
    if (el === null) {
      el = doc.createElement('style');
      el.setAttribute('data-u2-styles', '');
      doc.head.insertBefore(el, Styles._packageZone(doc.head));
    }
    else {
      const current = Number(el.getAttribute('data-revision'));
      if (current > revision || (current === revision && el.getAttribute('data-hash') === hash))
        return Styles.installed();
      if (current === revision) {
        if (!warned)
          console.warn(`u2: two copies of stylesheet revision ${revision} differ (${el.getAttribute('data-hash')} ` +
            `vs ${hash}) — a package was built from uncommitted u2 CSS; rebuild it.`);
        warned = true;
        return Styles.installed();
      }
    }
    el.textContent = css;
    el.setAttribute('data-revision', String(revision));
    el.setAttribute('data-hash', hash);
    return Styles.installed();
  }

  /** Where the sheet goes: after the platform's stylesheet links, before the first style a script
   * added after them — so the platform's rules come first, then u2's, then every package's, which
   * win a tie with u2's whatever order the packages load in. Null (the end) where there is none. */
  private static _packageZone(head: HTMLHeadElement): Element | null {
    const children = Array.from(head.children);
    const link = children.findIndex((el) => el.tagName === 'LINK' && el.getAttribute('rel') === 'stylesheet');
    return link < 0 ? null : children.slice(link + 1).find((el) => el.tagName === 'STYLE') ?? null;
  }

  /** The copies of u2's sheets an older build injected as plain styles: they come after the installed sheet and
   * win over it. Told by u2's own top-level `.u2-input-root` rule, which a package's overrides do not restate. */
  static legacyCopies(doc: Document = globalThis.document): number {
    return Array.from(doc.head.querySelectorAll('style')).filter((el) => !el.hasAttribute('data-u2-styles') &&
      /(^|}|\*\/)\s*\.u2-input-root\s*\{/.test(el.textContent ?? '')).length;
  }

  /** What the page runs, or null before any u2 copy installed. */
  static installed(): InstalledStyles | null {
    const el = globalThis.document?.head?.querySelector<HTMLStyleElement>(SELECTOR) ?? null;
    return el === null ? null :
      {revision: Number(el.getAttribute('data-revision')), hash: el.getAttribute('data-hash') ?? '', element: el};
  }
}

Styles.install();
