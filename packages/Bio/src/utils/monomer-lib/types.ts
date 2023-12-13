export class InvalidFilesError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'InvalidFilesError';
  }
}

