import sys
import json
import datetime
import websocket
from abc import ABC, abstractmethod

from .utils import DatagrokFanoutType
from .amqp_publisher import AmqpFanoutPublisher


class BaseRedirect(ABC):
    """Redirect class.

    Note:
        Redirects print statements. Should be only used with prefork workers and context manager (`with` statement).
    """

    service_log_key = "IS_SERVICE_LOG"
    max_message_count = 10000

    def __init__(self):
        self._original_stdout = sys.stdout
        self.message_count = 0

    @abstractmethod
    def write(self, message: str):
        pass 

    def formatMessage(self, message: str) -> dict:
        body = {
            "level": "debug",
            "time": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "message": message,
            "flag": "MISC",
            "params": {self.service_log_key: True}
        }     
        return body    

    def __enter__(self):
        sys.stdout = self
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout = self._original_stdout

    def flush(self):
        pass


class WebSocketRedirect(BaseRedirect):
    """Redirect stdout messages to WebSocket

    Note:
        Redirects print statements. Should be only used with prefork workers and context manager.
    """

    def __init__(self, ws: websocket.WebSocket):
        super().__init__()
        self.ws = ws
        self.notified_overflow = False

    def write(self, message: str):
        if self.message_count > self.max_message_count:
            if not self.notified_overflow:
                try:
                    self.ws.send(f"LOG {json.dumps({"level": "warning", 
                                                 "time": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                                                 "message": "Stdout overflow. Possible printing in the loop?",
                                                 "params": {self.service_log_key: True}})}")
                    self.notified_overflow = True
                except Exception:
                    pass
            return 
        try:
            self.message_count += 1
            self.ws.send(f"LOG {json.dumps(self.formatMessage(message))}") 
        except Exception:
            pass    
              

class AmqpRedirect(BaseRedirect):
    """Redirect stdout messages to Datagrok fanout

    Note:
        Redirects print statements. Should be only used with prefork workers and context manager.
    """

    def __init__(self, call_id: str):
        super().__init__()
        self.call_id = call_id
        self.amqp_publisher = AmqpFanoutPublisher.get_instance()
        self.notified_overflow = False
     
    def write(self, message: str):
        if self.message_count > self.max_message_count:
            if not self.notified_overflow:
                try:
                    self.amqp_publisher.publish({"level": "warning", 
                                                 "time": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                                                 "message": "Stdout overflow. Possible printing in the loop?",
                                                 "params": {self.service_log_key: True}}, self.call_id, DatagrokFanoutType.LOG)
                    self.notified_overflow = True
                except Exception:
                    pass
            return 
        try:   
            self.message_count += 1
            self.amqp_publisher.publish(self.formatMessage(message), self.call_id, DatagrokFanoutType.LOG)
        except Exception:
            pass                  
