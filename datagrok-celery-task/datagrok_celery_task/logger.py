import json
import logging
import re
import sys
import multiprocessing
from datetime import datetime, timezone

_logger = None


class SafeFormatter(logging.Formatter):
    def format(self, record):
        worker_name = multiprocessing.current_process().name
        record.worker_name = worker_name
        if not hasattr(record, 'task_id'):
            record.task_id = ''
        return super().format(record)


class JsonLogFormatter(logging.Formatter):
    """One JSON line per record in the shape the Datagrok server prints, plus `service`.
    Secrets are redacted as the server's `LogRedactor` does (Redaction contract, core/docs/AUDIT_LOGGING.md)."""
    REDACTED = '[REDACTED]'
    SECRET_KEYS = ('password', 'passwd', 'secret', 'token', 'apikey', 'authorization', 'cookie', 'credential',
                   'privatekey', 'accesskey')
    SECRETS = [re.compile(r'\beyJ[\w\-]+\.[\w\-]*\.[\w\-]*'), re.compile(r'\b(?:AKIA|ASIA)[A-Z0-9]{16}\b'),
               re.compile(r'-----BEGIN [A-Z ]*PRIVATE KEY-----[\s\S]*?(?:-----END [A-Z ]*PRIVATE KEY-----|\Z)')]
    PREFIXED = [re.compile(r'\b((?:Bearer|Basic)\s+)[\w\-.~+/]+=*', re.I),
                re.compile(r'(\b[a-zA-Z][\w+.\-]*://)[^\s/@:]+:[^\s/@]+(?=@)'),
                re.compile(r'''(\b[\w\-]*(?:password|passwd|pwd|secret|token|api[_\-]?key|access[_\-]?key)"?(?:\s*[=:]\s*|\s+(?=['"])))(?:"[^"]*"|'[^']*'|[^\s&,;"']+)''', re.I)]
    RECORD_ATTRS = set(vars(logging.makeLogRecord({}))) | {'message', 'asctime', 'taskName'}

    def __init__(self, service: str):
        super().__init__()
        self.service = service

    @classmethod
    def redact_text(cls, s):
        if not s:
            return s
        for r in cls.SECRETS:
            s = r.sub(cls.REDACTED, s)
        for r in cls.PREFIXED:
            s = r.sub(lambda m: m.group(1) + cls.REDACTED, s)
        return s

    @classmethod
    def redact_value(cls, v, key=None):
        if v is not None and key is not None and any(k in re.sub(r'[_\-]', '', str(key).lower()) for k in cls.SECRET_KEYS):
            return cls.REDACTED
        if isinstance(v, str):
            return cls.redact_text(v)
        if isinstance(v, dict):
            return {k: cls.redact_value(x, k) for k, x in v.items()}
        if isinstance(v, (list, tuple)):
            return [cls.redact_value(x) for x in v]
        return v

    def format(self, record):
        params = {k: v for k, v in vars(record).items() if k not in self.RECORD_ATTRS}
        line = {'v': 1,
                'time': datetime.fromtimestamp(record.created, timezone.utc).isoformat(timespec='milliseconds').replace('+00:00', 'Z'),
                'level': 'error' if record.levelno >= logging.ERROR else record.levelname.lower(),
                'service': self.service,
                'message': self.redact_text(record.getMessage())}
        if 'sessionId' in params:
            line['sessionId'] = params.pop('sessionId')
        line['params'] = self.redact_value(params)
        if record.exc_info:
            line['stackTrace'] = self.redact_text(self.formatException(record.exc_info))
        return json.dumps(line, default=str)


def setup_logger(name: str, level: int = logging.INFO, json_logs: bool = False) -> logging.Logger:
    logger = logging.getLogger(name)
    logger.setLevel(level)

    if not logger.handlers:
        handler = logging.StreamHandler(sys.stdout)
        handler.setLevel(level)
        handler.setFormatter(JsonLogFormatter(name) if json_logs else
                             SafeFormatter("[%(asctime)s: %(levelname)s/%(worker_name)s] %(task_id)s: %(message)s"))
        logger.addHandler(handler)
        logger.propagate = False

    return logger


def get_logger():
    global _logger
    if _logger is None:
        from .settings import Settings
        settings = Settings.get_instance()
        _logger = setup_logger(settings.celery_name, level=settings.log_level, json_logs=settings.json_logs)
    return _logger
