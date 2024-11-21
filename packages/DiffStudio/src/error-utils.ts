// Tools for Diff Studio model errors

export class ModelError extends Error {
  public helpUrl: string;

  constructor(message: string, helpUrl: string) {
    super(message);
    this.helpUrl = helpUrl;
  }
};
