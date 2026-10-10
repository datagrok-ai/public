/** Minimal timestamped stdout logger (mirrors datagrok-celery-task's logger.py role). */

const REDACTED = '[REDACTED]';
const SECRETS = [/\beyJ[\w-]+\.[\w-]*\.[\w-]*/g, /\b(?:AKIA|ASIA)[A-Z0-9]{16}\b/g,
  /-----BEGIN [A-Z ]*PRIVATE KEY-----[\s\S]*?(?:-----END [A-Z ]*PRIVATE KEY-----|$)/g];
const PREFIXED = [/\b((?:Bearer|Basic)\s+)[\w\-.~+/]+=*/gi, /(\b[a-zA-Z][\w+.\-]*:\/\/)[^\s/@:]+:[^\s/@]+(?=@)/g,
  /(\b[\w\-]*(?:password|passwd|pwd|secret|token|api[_\-]?key|access[_\-]?key)"?(?:\s*[=:]\s*|\s+(?=['"])))(?:"[^"]*"|'[^']*'|[^\s&,;"']+)/gi];

/** Redacts secrets as the Datagrok server's `LogRedactor` does (Redaction contract, core/docs/AUDIT_LOGGING.md). */
export function redactText(s: string): string {
  for (const re of SECRETS)
    s = s.replace(re, REDACTED);
  for (const re of PREFIXED)
    s = s.replace(re, (_, kept) => kept + REDACTED);
  return s;
}

/** With GROK_JSON_LOGS=true, one JSON line in the shape the Datagrok server prints, plus `service`. */
export function line(level: string, message: string, taskId?: string, env: NodeJS.ProcessEnv = process.env): string {
  if ((env['GROK_JSON_LOGS'] ?? '').toLowerCase() === 'true') {
    return JSON.stringify({v: 1, time: new Date().toISOString(), level: level === 'WARN' ? 'warning' : level.toLowerCase(),
      service: env['DATAGROK_CELERY_NAME'] || 'datagrok-celery', message: redactText(message),
      params: taskId ? {taskId} : {}});
  }
  return `${new Date().toISOString()} [${level}]${taskId ? ` [task ${taskId}]` : ''} ${redactText(message)}`;
}

export function logInfo(message: string, taskId?: string): void {
  console.log(line('INFO', message, taskId));
}

export function logWarn(message: string, taskId?: string): void {
  console.log(line('WARN', message, taskId));
}

export function logError(message: string, taskId?: string): void {
  console.error(line('ERROR', message, taskId));
}
