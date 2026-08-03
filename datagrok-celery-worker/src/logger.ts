/** Minimal timestamped stdout logger (mirrors datagrok-celery-task's logger.py role). */

function line(level: string, message: string, taskId?: string): string {
  return `${new Date().toISOString()} [${level}]${taskId ? ` [task ${taskId}]` : ''} ${message}`;
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
