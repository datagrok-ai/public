
export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}
