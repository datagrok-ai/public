import time as time_lib
from typing import Any, Dict, Optional
import uuid
import json
import inspect
import traceback
import threading
import signal
import sys

import celery
import websocket
from celery.signals import worker_shutdown, task_revoked
from celery import current_app

from .utils import InputValueProcessor, ReturnValueProcessor, camel_to_snake, DatagrokFanoutType
from .func_call import FuncCall, Type, FuncCallParam, FuncCallStatus
from .settings import Settings
from .logger import get_logger
from .redirect import WebSocketRedirect, AmqpRedirect
from .amqp_publisher import AmqpFanoutPublisher


class FuncCallRegistry:
    """Registry of currently tracked FuncCall objects keyed by their unique task ID. 
    Thread safe to support other than prefork pools.
    In a prefork worker pool it won't work since we terminate the whole process.
    """
    def __init__(self):
        self._lock = threading.Lock()
        self._registry: Dict[str, FuncCall] = {}

    def add(self, task_id: str, func_call: FuncCall) -> None:
        with self._lock:
            self._registry[task_id] = func_call

    def get(self, task_id: str) -> Optional[FuncCall]:
        with self._lock:
            return self._registry.get(task_id)

    def remove(self, task_id: str) -> None:
        with self._lock:
            self._registry.pop(task_id, None)

    def contains(self, task_id: str) -> bool:
        with self._lock:
            return task_id in self._registry
        
        
_task_registry = FuncCallRegistry()        


@task_revoked.connect
def on_task_cancel(request, terminated, signum, expired, **kwargs):
    task_id = getattr(request, 'id', None)
    if task_id is None:
        return

    if not _task_registry.contains(task_id):
        return

    call = _task_registry.get(task_id)
    send_canceled_call(call)  


@worker_shutdown.connect
def on_worker_shutdown(**kwargs):
    if AmqpFanoutPublisher._initialized:
        try:
            AmqpFanoutPublisher.get_instance().close()
        except Exception:
            pass 


def get_worker_pool_type() -> str:
    return current_app.conf.worker_pool


def send_canceled_call(call: FuncCall):   
    if call.status in (FuncCallStatus.COMPLETED, FuncCallStatus.ERROR, FuncCallStatus.CANCELED):
        return    
    try:
        call.status = FuncCallStatus.CANCELED
        amqp_publisher = AmqpFanoutPublisher.get_instance()
        amqp_publisher.publish(call.to_json(), call.id, DatagrokFanoutType.CALL)
        _task_registry.remove(call.id)
    except Exception:
        pass  


class DatagrokTask(celery.Task):
    # Disable Celery argument checking
    typing = False

    def __init__(self):
        self._settings = Settings.get_instance()
        self._logger = get_logger()
        self._pipe = None
        self._is_prefork = get_worker_pool_type() == 'prefork'

    # Overrides celery.Task.before_start
    def before_start(self, task_id, args, kwargs):
        """This function can't throw. It is special hook called before the tasks lifecycle starts, so exceptions here won't be catched 
        and task will never start.
        """
        if AmqpFanoutPublisher._initialized:
            self._amqp_publisher = AmqpFanoutPublisher.get_instance()
        else:
            self._amqp_publisher = AmqpFanoutPublisher(self._settings.broker_url, self._settings.calls_fanout)

        # Register revoke signal for prefork worker
        if self._is_prefork:
            def _on_cancel(signum, frame):
                if self._call is not None:
                    send_canceled_call(self._call)
                sys.exit(0)
            signal.signal(signal.SIGTERM, _on_cancel)

        self._call: FuncCall = FuncCall(args[0])
        self._accepted_send = False 
        _task_registry.add(self._call.id, self._call)
        # Try to establish AMQP connection to notify that call was accepted
        try:
            self._amqp_publisher.publish({}, self._call.id, DatagrokFanoutType.ACCEPTED)
            self._accepted_send = True
        except Exception as e:
            self._logger.warning("Initial Accepted message send failed: %s", str(e), extra={"task_id": self._call.id}, exc_info=True)   

        if self._call.requires_pipe:
            try:
                self._ensure_pipe_conn(timeout=30)     
            except Exception as e:
                self._pipe = None
                self._logger.warning("Initial grok_pipe connection failed with exception: %s", str(e), exc_info=True, extra={"task_id": self._call.id})


    # Overrides celery.Task.after_return
    def after_return(self, status, retval, task_id, args, kwargs, einfo):
        _task_registry.remove(self._call.id)
        # Unregister revoke signal for prefork worker
        if self._is_prefork:
            signal.signal(signal.SIGTERM, signal.SIG_DFL)

        self._call.status = FuncCallStatus.COMPLETED if status == 'SUCCESS' else FuncCallStatus.ERROR

        if einfo is not None:
            self._call.error_message = str(einfo.exception)
            self._call.error_stack_trace = "".join(einfo.traceback)
        elif isinstance(retval, Exception):
            self._call.error_message = str(retval)
            self._call.error_stack_trace = ''.join(traceback.format_exception(type(retval), retval, retval.__traceback__))

        for i in range(3):
            if i > 0:
                self._logger.info("Retrying (%d) CALL send", i + 1, extra={"task_id": self._call.id})
            try:
                if self._call.requires_pipe:
                    self._pipe_send(f"CALL {json.dumps(self._call.to_json())}") 
                else:
                    self._amqp_publisher.publish(self._call.to_json(), self._call.id, DatagrokFanoutType.CALL)
                break    
            except Exception as e:
                self._logger.error("CALL send failed (attempt %d): %s", i + 1, str(e), exc_info=True, extra={"task_id": self._call.id})
                time_lib.sleep(2 ** i) 

        if self._pipe is not None:
            try:
                self._pipe.close()
            except Exception as e:
                self._logger.error("Couldn't close WS connection: %s", str(e), extra={"task_id": self._call.id}, exc_info=True)
            finally:
                self._pipe = None
        self._call = None

    # Overrides celery.Task.update_state
    def update_state(self, task_id=None, state=None, meta=None):
        try:
            percent = meta.get("percent") if meta else None
            description = meta.get("description") if meta else None
            if percent is None:
                self._logger.warning("Incorrect usage of update_state, meta['percent'] key should be present.")
                return
            message = {'percent': percent, 'description': description}
            if self._call.requires_pipe:
                self._pipe_send(f"PROGRESS {json.dumps(message)}")
            else:
                self._amqp_publisher.publish(message, self._call.id, DatagrokFanoutType.PROGRESS)
        except Exception as e:
            # Ignoring update status errors
            self._logger.warning("Update status failed: %s", str(e), extra={"task_id": self._call.id}, exc_info=True)           

    # Overrides celery.Task.__call__
    def __call__(self, *args, **kwargs):
        if len(self._call.output_params) > 1:
            raise ValueError("Only one return parameter is allowed, but more are present in the Datagrok function declaration.") 
         
        if not self._accepted_send:
            try:
                self._amqp_publisher.publish({}, self._call.id, DatagrokFanoutType.ACCEPTED)
                self._accepted_send = True
            except Exception as e:
                self._logger.error("Accepted message send failed: %s", str(e), extra={"task_id": self._call.id}, exc_info=True)
                raise 
        if self._call.requires_pipe:
            self._ensure_pipe_conn(timeout=30)

        self._process_input_params()

        if not self._is_prefork:
            return self._run_and_process()

        if self._call.requires_pipe:
            with WebSocketRedirect(self._pipe):
                return self._run_and_process()
        else:
            with AmqpRedirect(self._call.id):
                return self._run_and_process()
   
    def _run_and_process(self):
        value = self._run()
        self._process_return_value(value)
        return None

    def _run(self) -> Any:
        args = []
        kwargs = {}
        param_map = {camel_to_snake(x.name): x for x in self._call.input_params}
        fn = getattr(self, '__wrapped__', None)
        if fn is None:
            raise RuntimeError("Incorrect usage of DatagrokTask. Should be used as a base for Celery decorators")
        
        sig = inspect.signature(fn)
        parameters = [x for x in sig.parameters.values()]
        prepend_self = (
            len(parameters) > 0 and
            parameters[0].name == 'self' and
            not inspect.ismethod(fn)  # method means it's already bound
        )
        if prepend_self:
            args.append(self)

        for param in parameters[prepend_self:]:
            name = param.name
            call_param = param_map.pop(name, None)
            if call_param is not None:
                if param.kind in (inspect.Parameter.POSITIONAL_ONLY, inspect.Parameter.POSITIONAL_OR_KEYWORD):
                    args.append(call_param.value)
                elif param.kind == inspect.Parameter.KEYWORD_ONLY:
                    kwargs[name] = call_param.value
            else:
                if param.kind == inspect.Parameter.VAR_POSITIONAL:
                    self._logger.warning("Celery task contains *args, skipping it", extra={"task_id": self._call.id})
                    continue
                elif param.kind == inspect.Parameter.VAR_KEYWORD:
                    # add session token that can be used in datagrok_api.DatagrokClient
                    param_map["USER_API_KEY"] = self._call.user_api_key
                    kwargs.update(param_map)    
                    break    
            
        try:
            sig.bind(*args, **kwargs)
        except TypeError as e:
            self._logger.error("Couldn't match Datagrok func parameters to the task: %s", str(e), extra={"task_id": self._call.id}, exc_info=True)
            raise

        return self.run(*args, **kwargs)
        
    def _process_return_value(self, value):
        outputs = self._call.output_params
        if not outputs:
            return

        param = outputs[0]

        if value is None:
            param.value = None
            self._logger.warning("Task returned None", extra={"task_id": self._call.id})
            return

        setter = ReturnValueProcessor.type_map.get(param.property_type)

        if setter is None:
            raise ValueError(f"Unsupported return type: {param.property_type}")
        
        if param.property_type == Type.DATA_FRAME:
            setter(param, value, self._call.use_parquet_transfer)
        else:
            setter(param, value)

        if param.is_streamable:
            self._send_param_grok_pipe(param)

    def _send_param_grok_pipe(self, param: FuncCallParam):
        self._logger.info("Sending param %s", param.name, extra={"task_id": self._call.id})
        if not isinstance(param.value, (bytes, bytearray)):
                raise TypeError(f"Expected bytes for param.value, got {type(param.value)}")
        
        param_type = 'blob'
        param_id = param.name
        data = param.value
        size = len(data)
        if param.property_type == Type.DATA_FRAME:
            param_type = "parquet" if self._call.use_parquet_transfer else "csv"
            param_id = str(uuid.uuid1())
        param.value = {"id": param_id} if param.property_type == Type.DATA_FRAME else param_id
        tags_json = json.dumps({".id": param_id, ".type": param_type})

        start = time_lib.monotonic()

        self._pipe_send(f"SENDING DATAFRAME {size} {tags_json}")

        batch_size = self._call.binary_batch_size
        for i in range(0, size, batch_size):
            chunk = data[i:i+batch_size]
            self._logger.debug("Sending binary batch %s of param %s", i, param.name, extra={"task_id": self._call.id})
            self._pipe_send(chunk, opcode=websocket.ABNF.OPCODE_BINARY)
            
            try:
                response = self._pipe.recv()
            except websocket.WebSocketTimeoutException as e:
                self._logger.error("Timeout waiting for message from grok_pipe: %s", str(e), exc_info=True, extra={"task_id": self._call.id})
                raise
            if response == "ERROR":
                raise Exception(f"Error sending parameter {param.name}")
            elif response != "PART OK":
                raise ValueError("Unexpected message received from grok_pipe.")
            
        self._logger.info("Finished sending param %s (%d bytes in %d chunks), elapsed: %d ms", param.name, size, (size + batch_size - 1) // batch_size, (time_lib.monotonic() - start) * 1000, extra={"task_id": self._call.id})

    def _process_input_params(self):
        for param in self._call.input_params:
            if param.is_streamable:
                param.value = self._get_param_grok_pipe(param)
            if param.value is not None:
                setter = InputValueProcessor.type_map.get(param.property_type)
                if setter is not None:
                    if param.property_type == Type.DATA_FRAME:
                        setter(param, self._call.use_parquet_transfer)
                    else:
                        setter(param)
                  
    def _get_param_grok_pipe(self, param: FuncCallParam):
        self._logger.info("Receiving param %s", param.name, extra={"task_id": self._call.id})

        start = time_lib.monotonic()

        self._pipe_send(f"PARAM {param.name}")
        
        chunks = []
        total_bytes_received = 0
        expected_size = None
        
        while True:
            if (time_lib.monotonic() - start) > self._settings.param_timeout_minutes * 60:
                self._logger.error("Timeout receiving param %s", param.name, extra={"task_id": self._call.id})
                raise TimeoutError(f"Timeout receiving param {param.name}")

            try:
                message = self._pipe.recv()
            except websocket.WebSocketTimeoutException as e:
                self._logger.error("Timeout waiting for message from grok_pipe: %s", str(e), exc_info=True, extra={"task_id": self._call.id})
                raise
                
            if isinstance(message, str):
                if message.startswith(f"PARAM_SENT {param.name}"):
                    break
                elif message.startswith("SENDING"):
                    self._logger.debug("Server is sending parameter data: %s", message, extra={"task_id": self._call.id})
                    expected_size = self._get_expected_size(message)
                    self._logger.debug("Expected parameter size: %s bytes", expected_size)
                else:
                    self._logger.debug(f"Unexpected message: %s", message[:30])
            else:
                if expected_size is None:
                    self._logger.error("Expected size was not set before binary data began for param %s.", param.name, extra={"task_id": self._call.id})
                    raise ValueError("Protocol error: expected size not defined.")
                
                chunk_size = len(message)
                self._logger.debug("Received binary chunk of size %s", chunk_size, extra={"task_id": self._call.id})
                chunks.append(message)
                total_bytes_received += chunk_size
                self._logger.debug("Total bytes received: %s/%s", total_bytes_received, expected_size, extra={"task_id": self._call.id})
                self._pipe_send("PART OK")
                if total_bytes_received > expected_size:
                    self._logger.error("Received size of binary data is larger than expected for param %s", param.name, extra={"task_id": self._call.id})
                    raise ValueError("Received size of binary data is larger than expected.")      
                if total_bytes_received == expected_size:
                    self._logger.info("Received all bytes for param %s, elapsed: %d ms", param.name, (time_lib.monotonic() - start) * 1000, extra={"task_id": self._call.id})
        
        if not chunks:
            return None
        elif len(chunks) == 1:
            return chunks[0]
        else:
            return b''.join(chunks)

    def _get_expected_size(self, message: str) -> int:
        try:
            parts = message.split()
            if len(parts) >= 3:
                return int(parts[2])
        except (ValueError, IndexError) as e:
            self._logger.error("Could not parse size from SENDING: %s", str(e), exc_info=True, extra={"task_id": self._call.id})
            raise

    def _ensure_pipe_conn(self, max_retries=5, retry_interval=3, timeout=15):
        if self._pipe is not None and self._pipe.connected:
            return
        attempt = 0
        while attempt < max_retries:
            try:
                self._pipe = websocket.create_connection(f"{self._settings.pipe_url}/{self._call.id}", header=[f"x-member-name: celery-{self._settings.celery_name}", f"authorization: {self._settings.pipe_key}"] or [], timeout=timeout)
                self._pipe.settimeout(self._settings.ws_message_timeout_seconds)
                return
            except Exception as e:
                self._logger.warning("WS connection failed (attempt %s/%s): %s", attempt + 1, max_retries, str(e), exc_info=True, extra={"task_id": self._call.id})
                attempt += 1
                time_lib.sleep(retry_interval)
        raise ConnectionError(f"Failed to connect to grok_pipe at {self._settings.pipe_url} after {max_retries} retries") 

    def _pipe_send(self, message, opcode=websocket.ABNF.OPCODE_TEXT):
        self._ensure_pipe_conn()
        try:
            self._pipe.send(message, opcode=opcode)
        except websocket.WebSocketConnectionClosedException:
            self._logger.warning("Grok_Pipe is disconnected, trying to reconnect...", self._call.id, extra={"task_id": self._call.id})
            self._ensure_pipe_conn()
            self._pipe.send(message) 
            