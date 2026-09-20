/* `?raw` imports (rspack `asset/source`, see rspack.config.js) resolve to the file's text. */
declare module '*?raw' {
  const text: string;
  export default text;
}
