export const WORKSPACE = process.env['CLAUDE_WORKSPACE'] || '/workspace';

export function rewriteForDocker(url: string): string {
  return url.replace('//localhost', '//host.docker.internal');
}
