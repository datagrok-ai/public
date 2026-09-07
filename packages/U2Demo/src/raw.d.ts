/* `?raw` imports (webpack `asset/source`, see webpack.config.js) resolve to the file's text. */
declare module '*?raw' {
  const text: string;
  export default text;
}
