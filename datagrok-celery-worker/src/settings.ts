/** Worker settings, read from environment variables.
 *  Mirrors datagrok-celery-task/settings.py (including the existing 'AMPQ' spelling
 *  of the AMQP variables — kept for parity with the python lib and deploy tooling). */
export class Settings {
  taskQueueName: string;
  celeryHostname: string;
  celeryName: string;
  packageName: string;
  packageVersion?: string;
  amqpHost: string;
  amqpPort: number;
  amqpUser: string;
  amqpPassword: string;
  amqpTls: boolean;
  pipeHost: string;
  pipePort: number;
  pipeKey: string;
  apiUrl?: string;
  paramTimeoutMinutes: number;
  wsMessageTimeoutSeconds: number;
  healthPort: number;

  constructor(values: {
    taskQueueName: string,
    celeryHostname: string,
    celeryName: string,
    packageName: string,
    packageVersion?: string,
    amqpHost: string,
    amqpPort: number,
    amqpUser: string,
    amqpPassword: string,
    amqpTls: boolean,
    pipeHost: string,
    pipePort: number,
    pipeKey: string,
    apiUrl?: string,
    paramTimeoutMinutes: number,
    wsMessageTimeoutSeconds: number,
    healthPort: number,
  }) {
    this.taskQueueName = values.taskQueueName;
    this.celeryHostname = values.celeryHostname;
    this.celeryName = values.celeryName;
    this.packageName = values.packageName;
    this.packageVersion = values.packageVersion;
    this.amqpHost = values.amqpHost;
    this.amqpPort = values.amqpPort;
    this.amqpUser = values.amqpUser;
    this.amqpPassword = values.amqpPassword;
    this.amqpTls = values.amqpTls;
    this.pipeHost = values.pipeHost;
    this.pipePort = values.pipePort;
    this.pipeKey = values.pipeKey;
    this.apiUrl = values.apiUrl;
    this.paramTimeoutMinutes = values.paramTimeoutMinutes;
    this.wsMessageTimeoutSeconds = values.wsMessageTimeoutSeconds;
    this.healthPort = values.healthPort;
  }

  static fromEnv(env: NodeJS.ProcessEnv = process.env): Settings {
    const missing: string[] = [];
    const required = (name: string): string => {
      const value = env[name];
      if (value == null || value === '')
        missing.push(name);
      return value ?? '';
    };
    const intEnv = (name: string, defaultValue: number): number => {
      const parsed = parseInt(env[name] ?? '', 10);
      return Number.isFinite(parsed) ? parsed : defaultValue;
    };

    const settings = new Settings({
      taskQueueName: required('TASK_QUEUE_NAME'),
      celeryHostname: required('CELERY_HOSTNAME'),
      celeryName: env['DATAGROK_CELERY_NAME'] || 'datagrok-celery',
      packageName: required('DATAGROK_PACKAGE_NAME'),
      packageVersion: env['DATAGROK_PACKAGE_VERSION'] || undefined,
      amqpHost: env['DATAGROK_AMPQ_HOST'] || 'localhost',
      amqpPort: intEnv('DATAGROK_AMPQ_PORT', 5672),
      amqpUser: env['DATAGROK_AMPQ_USER'] || 'guest',
      amqpPassword: env['DATAGROK_AMPQ_PASSWORD'] || 'guest',
      amqpTls: (env['DATAGROK_AMPQ_TLS'] ?? 'false').toLowerCase() === 'true',
      pipeHost: env['DATAGROK_PIPE_HOST'] || 'localhost',
      pipePort: intEnv('DATAGROK_PIPE_PORT', 3000),
      pipeKey: env['DATAGROK_PIPE_KEY'] ?? '',
      apiUrl: env['DATAGROK_API_URL'] || undefined,
      paramTimeoutMinutes: intEnv('DATAGROK_PARAM_TIMEOUT', 5),
      wsMessageTimeoutSeconds: intEnv('DATAGROK_WS_MESSAGE_TIMEOUT', 30),
      healthPort: intEnv('HEALTH_PORT', 8000),
    });
    if (missing.length > 0)
      throw new Error(`Missing required environment variables: ${missing.join(', ')}`);
    return settings;
  }

  /** amqp(s)://user:password@host:port (settings.py builds the same URL; credentials
   *  are URL-encoded here because amqplib parses the string as a URL). */
  get brokerUrl(): string {
    const protocol = this.amqpTls ? 'amqps' : 'amqp';
    return `${protocol}://${encodeURIComponent(this.amqpUser)}:${encodeURIComponent(this.amqpPassword)}@${this.amqpHost}:${this.amqpPort}`;
  }

  get pipeUrl(): string {
    return `ws://${this.pipeHost}:${this.pipePort}`;
  }
}
