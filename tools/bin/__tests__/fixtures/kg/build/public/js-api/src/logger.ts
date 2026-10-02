export class Logger {
  static getStatic(): Logger {
    return new Logger();
  }

  info(message: string): void {
  }
}
