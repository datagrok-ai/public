import logging
import sys
import multiprocessing

_logger = None


class SafeFormatter(logging.Formatter):
    def format(self, record):
        worker_name = multiprocessing.current_process().name
        record.worker_name = worker_name
        if not hasattr(record, 'task_id'):
            record.task_id = ''
        return super().format(record)


def setup_logger(name: str, level: int = logging.INFO) -> logging.Logger:
    logger = logging.getLogger(name)
    logger.setLevel(level)

    if not logger.handlers:
        handler = logging.StreamHandler(sys.stdout)
        handler.setLevel(level)
        handler.setFormatter(SafeFormatter("[%(asctime)s: %(levelname)s/%(worker_name)s] %(task_id)s: %(message)s"))
        logger.addHandler(handler)
        logger.propagate = False

    return logger


def get_logger():
    global _logger
    if _logger is None:
        from .settings import Settings
        settings = Settings.get_instance()
        _logger = setup_logger(settings.celery_name, level=settings.log_level)
    return _logger
