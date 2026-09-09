import {CommandMessage, PatinaeViewer} from './patinae-types';

export type PmlLine =
  {kind: 'include', path: string} |
  {kind: 'load', path: string, name?: string} |
  {kind: 'command', text: string};

// Mirrors the upstream script engine (crates/patinae-cmd/src/script.rs): backslash continuations,
// blank and # lines skipped, @file includes, load lines with an optional object name.
export function parsePml(text: string): PmlLine[] {
  const lines = text.replace(/\s*\\\r?\n\s*/g, ' ').split(/\r?\n/);
  const result: PmlLine[] = [];
  for (const raw of lines) {
    const line = raw.trim();
    if (line === '' || line.startsWith('#'))
      continue;
    if (line.startsWith('@')) {
      result.push({kind: 'include', path: line.slice(1).trim()});
      continue;
    }
    const load = line.match(/^load\s+([^,]+?)\s*(?:,\s*([^,]+?)\s*)?(?:,.*)?$/i);
    if (load && !/^https?:\/\//i.test(load[1]))
      result.push({kind: 'load', path: load[1], name: load[2]});
    else
      result.push({kind: 'command', text: line});
  }
  return result;
}

export interface PmlFileReader {
  readText(path: string): Promise<string>;
  readBytes(path: string): Promise<Uint8Array>;
}

export function baseName(path: string): string {
  return path.split('/').pop()!.replace(/\.gz$/i, '').replace(/\.[^.]+$/, '');
}

export function extension(path: string): string {
  return path.replace(/\.gz$/i, '').split('.').pop()!.toLowerCase();
}

export type PmlEcho = (command: string, messages: CommandMessage[]) => void;

// Runs a script whose relative paths resolve against `dir` (a Datagrok file-share folder);
// `echo` receives every executed command with its output, e.g. to transcribe it into the command panel.
export async function runPml(text: string, dir: string, files: PmlFileReader, viewer: PatinaeViewer, echo?: PmlEcho):
  Promise<CommandMessage[]> {
  const messages: CommandMessage[] = [];
  for (const line of parsePml(text)) {
    if (line.kind === 'include')
      messages.push(...await runPml(await files.readText(`${dir}/${line.path}`), dir, files, viewer, echo));
    else if (line.kind === 'load') {
      const name = line.name ?? baseName(line.path);
      viewer.loadData(await files.readBytes(`${dir}/${line.path}`), name, extension(line.path));
      echo?.(`load ${line.path}, ${name}`, [{level: 'info', text: ` Loaded "${name}"`}]);
    } else {
      const output = (await viewer.executeAsync(line.text)).messages;
      messages.push(...output);
      echo?.(line.text, output);
    }
  }
  return messages;
}
