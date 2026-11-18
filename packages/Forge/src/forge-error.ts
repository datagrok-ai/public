export class ForgeError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'ForgeError';
  }
}

export function errorMessage(e: unknown): string {
  return e instanceof Error ? e.message : String(e);
}
