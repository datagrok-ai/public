import time
from threading import Lock
from typing import Optional, Dict, Any
from kombu import Connection, Exchange, Producer
from amqp.exceptions import (
    ConnectionError,
    ChannelError
)
from kombu.exceptions import KombuError

from .utils import DatagrokFanoutType


class AmqpFanoutPublisher:
    _instance: Optional["AmqpFanoutPublisher"] = None
    _lock = Lock()

    _initialized: bool = False

    _broker_url: str
    _exchange_name: str
    _durable: bool

    _amqp_conn: Optional[Connection] = None
    _exchange: Optional[Exchange] = None
    _producer: Optional[Producer] = None

    _max_retries: int = 3
    _connection_timeout: float = 10.0
    _publish_timeout: float = 10.0

    def __new__(cls, *args, **kwargs):
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self, broker_url: str, exchange_name: str, durable: bool = False):
        if AmqpFanoutPublisher._initialized:
            return

        self._broker_url = broker_url
        self._exchange_name = exchange_name
        self._durable = durable

        AmqpFanoutPublisher._initialized = True

    @classmethod
    def get_instance(cls) -> "AmqpFanoutPublisher":
        if cls._instance is None or not cls._initialized:
            raise RuntimeError("AmqpFanoutPublisher is not initialized yet. "
                               "Create an instance first with broker_url and exchange_name.")
        return cls._instance

    def _ensure_amqp_conn(self) -> None:
        retry_count = 0
        last_error = None

        while retry_count < self._max_retries:
            try:
                if self._amqp_conn:
                    try:
                        self._amqp_conn.ensure_connection(max_retries=1, timeout=self._connection_timeout)
                        return
                    except KombuError:
                        self._amqp_conn = None

                self._amqp_conn = Connection(self._broker_url, heartbeat=60)
                self._amqp_conn.connect()

                self._exchange = Exchange(
                    self._exchange_name,
                    type="fanout",
                    durable=self._durable
                )
                self._exchange(self._amqp_conn).declare()
                self._producer = Producer(self._amqp_conn)
                return

            except KombuError as e:
                last_error = e
                retry_count += 1
                if retry_count < self._max_retries:
                    time.sleep(2 ** retry_count)

                self._amqp_conn = None
                self._exchange = None
                self._producer = None

        raise KombuError(f"Failed to establish AMQP connection after "
                         f"{self._max_retries} attempts: {last_error}")

    def publish(
        self,
        payload: Dict[str, Any],
        correlation_id: str,
        type_: DatagrokFanoutType,
        routing_key: str = '',
        serializer: str = 'json'
    ) -> None:
        if not AmqpFanoutPublisher._initialized:
            raise RuntimeError("AmqpFanoutPublisher is not initialized")

        try:
            self._ensure_amqp_conn()
            self._producer.publish(
                payload,
                exchange=self._exchange,
                routing_key=routing_key,
                correlation_id=correlation_id,
                type=type_.value,
                serializer=serializer,
                timeout=self._publish_timeout
            )
        except (ConnectionError, ChannelError):
            self._amqp_conn = None
            self._producer = None
            self._exchange = None
            raise


    def close(self) -> None:
        if self._amqp_conn:
            try:
                self._amqp_conn.release()
            except Exception:
                pass
            finally:
                self._amqp_conn = None
                self._exchange = None
                self._producer = None
                AmqpFanoutPublisher._initialized = False
                AmqpFanoutPublisher._instance = None
