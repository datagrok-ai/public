
export interface ILogger {
  error(message: any, params?: object | undefined, stackTrace?: string | undefined): void;
  warning(message: string, params?: object | undefined): void;
  info(message: string, params?: object | undefined): void;
  debug(message: string, params?: object | undefined): void;
}
