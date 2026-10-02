export interface ParsedId {
  scheme: string;
  local: string;
}

export class Ids {
  parse(id: string): string {
    return id.split(':')[0];
  }
}

export function fileId(file: string): string {
  return `file:${file}`;
}
