/* `?raw` imports (rspack `asset/source`, see rspack.config.js) resolve to the file's text. */
declare module '*?raw' {
  const text: string;
  export default text;
}

/* rspack's `import.meta.webpackContext` (src/source-panel.ts bundles the page sources through one). */
interface ImportMeta {
  webpackContext(request: string, options?: {recursive?: boolean, regExp?: RegExp}): {
    keys(): string[];
    <T>(key: string): T;
  };
}
