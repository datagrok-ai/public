/* eslint-disable */
/**
 * Polyfills for browser APIs that this plugin uses but that are missing in Chrome ≤ 53 (the engine
 * used by Dartium, which is a regular test target for Grokky). Imported as the first import from
 * `package.ts` / `package-test.ts`, so these are in place before any other module runs.
 *
 * Keep this minimal — add an entry only when the plugin actually starts using the API.
 * The bundle targets `es6`, and Chrome 50 supports ES2015 syntax; newer syntax (`?.`, `??`,
 * object spread, `async`/`await`) is transpiled by TypeScript, so syntax is not a concern here —
 * only runtime APIs are.
 */

(function installPolyfills() {
  // crypto.randomUUID — Chrome 92+
  const cryptoObj: any = (typeof crypto !== 'undefined') ? crypto : undefined;
  if (cryptoObj && typeof cryptoObj.randomUUID !== 'function') {
    try {
      cryptoObj.randomUUID = function randomUUID(): string {
        return 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, (c) => {
          const r = Math.random() * 16 | 0;
          return (c === 'x' ? r : (r & 0x3 | 0x8)).toString(16);
        });
      };
    } catch (e) {
      console.warn('Grokky: could not polyfill crypto.randomUUID', e);
    }
  }

  // Promise.prototype.finally — Chrome 63+
  const P: any = Promise.prototype;
  if (!P.finally)
    P.finally = function (this: Promise<any>, cb: () => void) {
      return this.then(
        (v) => Promise.resolve(cb()).then(() => v),
        (e) => Promise.resolve(cb()).then(() => { throw e; }),
      );
    };

  // Object.values / Object.entries — Chrome 54+
  const O: any = Object;
  if (!O.values)
    O.values = (o: any) => Object.keys(o).map((k) => o[k]);
  if (!O.entries)
    O.entries = (o: any) => Object.keys(o).map((k) => [k, o[k]]);
  // Object.fromEntries — Chrome 73+
  if (!O.fromEntries)
    O.fromEntries = (it: Iterable<[any, any]>) => {
      const r: any = {};
      for (const [k, v] of it)
        r[k] = v;
      return r;
    };

  // String.prototype.trimStart / trimEnd — Chrome 66+ (Chrome 50 has the non-standard trimLeft/trimRight)
  const S: any = String.prototype;
  if (!S.trimStart)
    S.trimStart = S.trimLeft || function (this: string) { return this.replace(/^\s+/, ''); };
  if (!S.trimEnd)
    S.trimEnd = S.trimRight || function (this: string) { return this.replace(/\s+$/, ''); };

  // String.prototype.matchAll — Chrome 73+ (returns an array, which for..of iterates the same way)
  if (!S.matchAll)
    S.matchAll = function (this: string, re: RegExp) {
      const g = re.global ? re : new RegExp(re.source, re.flags + 'g');
      const out: RegExpExecArray[] = [];
      let m: RegExpExecArray | null;
      while ((m = g.exec(this)) !== null) {
        out.push(m);
        if (m[0] === '')
          g.lastIndex++;
      }
      return out;
    };

  // Array.prototype.flatMap — Chrome 69+ (Array.prototype.concat already flattens one level of array args)
  const A: any = Array.prototype;
  if (!A.flatMap)
    A.flatMap = function (this: any[], fn: (v: any, i: number, a: any[]) => any, thisArg?: any) {
      return this.reduce((acc: any[], v, i, a) => acc.concat(fn.call(thisArg, v, i, a)), []);
    };

  // ChildNode.replaceWith / ParentNode.append / ParentNode.prepend — Chrome 54+
  const E: any = (typeof Element !== 'undefined') ? Element.prototype : undefined;
  if (E) {
    const toNode = (x: any): Node => x instanceof Node ? x : document.createTextNode(String(x));
    if (!E.append)
      E.append = function (this: Element, ...args: any[]) {
        for (const a of args) this.appendChild(toNode(a));
      };
    if (!E.prepend)
      E.prepend = function (this: Element, ...args: any[]) {
        const ref = this.firstChild;
        for (const a of args) this.insertBefore(toNode(a), ref);
      };
    if (!E.replaceWith)
      E.replaceWith = function (this: Element, ...args: any[]) {
        if (!this.parentNode) return;
        const frag = document.createDocumentFragment();
        for (const a of args) frag.appendChild(toNode(a));
        this.parentNode.replaceChild(frag, this);
      };
  }
})();
