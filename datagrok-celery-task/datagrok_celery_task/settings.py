import logging
import os
from typing import Optional, ClassVar


class Settings:
    _instance: ClassVar[Optional["Settings"]] = None
    _initialized: ClassVar[bool] = False

    celery_name: str
    amqp_host: str
    amqp_user: str
    amqp_password: str
    amqp_port: int
    tls: bool
    pipe_host: str
    pipe_key: Optional[str]
    pipe_port: int
    param_timeout_minutes: int
    ws_message_timeout_seconds: int
    calls_fanout: str
    log_level: int
    broker_url: str
    pipe_url: str

    def __new__(cls, *args, **kwargs):
        if cls._instance is None:
            cls._instance = super(Settings, cls).__new__(cls)
        return cls._instance

    def __init__(
        self,
        celery_name: Optional[str] = None,
        amqp_host: Optional[str] = None,
        amqp_user: Optional[str] = None,
        amqp_password: Optional[str] = None,
        amqp_port: Optional[int] = None,
        amqp_tls: Optional[bool] = None,
        pipe_host: Optional[str] = None,
        pipe_key: Optional[str] = None,
        pipe_port: Optional[int] = None,
        param_timeout_minutes: Optional[int] = None,
        ws_message_timeout_seconds: Optional[int] = None,
        calls_fanout: Optional[str] = "calls_fanout",
        log_level: Optional[int] = logging.INFO
    ):
        if self._initialized:
            return

        self.celery_name = celery_name or os.getenv("DATAGROK_CELERY_NAME", "datagrok-celery")
        self.amqp_host = amqp_host or os.getenv("DATAGROK_AMPQ_HOST", "localhost")
        self.amqp_user = amqp_user or os.getenv("DATAGROK_AMPQ_USER", "guest")
        self.amqp_password = amqp_password or os.getenv("DATAGROK_AMPQ_PASSWORD", "guest")
        self.amqp_port = amqp_port if amqp_port is not None else self._get_int_env("DATAGROK_AMPQ_PORT", 5672)
        self.tls = amqp_tls if amqp_tls is not None else self._get_bool_env("DATAGROK_AMPQ_TLS", False)

        self.pipe_host = pipe_host or os.getenv("DATAGROK_PIPE_HOST", "localhost")
        self.pipe_key = pipe_key or os.getenv("DATAGROK_PIPE_KEY")
        self.pipe_port = pipe_port if pipe_port is not None else self._get_int_env("DATAGROK_PIPE_PORT", 3000)

        self.param_timeout_minutes = param_timeout_minutes if param_timeout_minutes is not None else self._get_int_env("DATAGROK_PARAM_TIMEOUT", 5)
        self.ws_message_timeout_seconds = ws_message_timeout_seconds if ws_message_timeout_seconds is not None else self._get_int_env("DATAGROK_WS_MESSAGE_TIMEOUT", 30)

        self.calls_fanout = calls_fanout
        self.log_level = log_level

        self.broker_url = self._build_broker_url()
        self.pipe_url = self._build_pipe_url()

        self._validate_settings()
        
        self._initialized = True

    @classmethod
    def get_instance(cls) -> "Settings":
        if cls._instance is None:
            raise RuntimeError("Settings must be initialized before getting instance")
        return cls._instance

    @staticmethod
    def _get_int_env(key: str, default: int) -> int:
        try:
            return int(os.getenv(key, str(default)))
        except (ValueError, TypeError):
      
            return default

    @staticmethod
    def _get_bool_env(key: str, default: bool = False) -> bool:
        return os.getenv(key, str(default)).lower() == "true"

    def _build_broker_url(self) -> str:
        protocol = "amqps" if self.tls else "amqp"
        return f"{protocol}://{self.amqp_user}:{self.amqp_password}@{self.amqp_host}:{self.amqp_port}"

    def _build_pipe_url(self) -> str:
        return f"ws://{self.pipe_host}:{self.pipe_port}"

    def _validate_settings(self) -> None:
        if not self.celery_name:
            raise ValueError("celery_name cannot be empty")
        if not self.amqp_host:
            raise ValueError("amqp_host cannot be empty")
        if not self.amqp_user:
            raise ValueError("amqp_user cannot be empty")
        if not self.amqp_password:
            raise ValueError("amqp_password cannot be empty")
        if self.amqp_port <= 0:
            raise ValueError("amqp_port must be positive")
        if not self.pipe_host:
            raise ValueError("pipe_host cannot be empty")
        if self.pipe_port <= 0:
            raise ValueError("pipe_port must be positive")
        if self.param_timeout_minutes <= 0:
            raise ValueError("param_timeout_minutes must be positive")
        if self.ws_message_timeout_seconds <= 0:
            raise ValueError("ws_message_timeout_seconds must be positive")
    