declare module 'aws4' {
  interface SignOptions {
    host?: string;
    method?: string;
    path?: string;
    service?: string;
    region?: string;
    headers?: Record<string, string | string[] | undefined>;
    body?: string | Buffer;
  }
  interface Credentials {
    accessKeyId: string;
    secretAccessKey: string;
    sessionToken?: string;
  }
  export function sign(options: SignOptions, credentials?: Credentials): SignOptions;
}
