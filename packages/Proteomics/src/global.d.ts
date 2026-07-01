// Ambient module declarations for webpack `asset/source` raw-text imports.
// Used by src/tests/spectronaut-parser.ts to inline the committed
// tools/spectronaut-aggregate.{sql,sh} as strings (Plan 12-03 R5 drift guard);
// the tools/ dir is not deployed under files/, so it is unreachable via
// _package.files at runtime — webpack asset/source is the only way to pin the
// committed fallback content into the test bundle.
declare module '*.sql' {
  const content: string;
  export default content;
}

declare module '*.sh' {
  const content: string;
  export default content;
}
