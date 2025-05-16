import time as time_lib
from typing import Any
import uuid
import json
import inspect
import traceback

import celery
import websocket
from celery.signals import worker_shutdown

from .utils import InputValueProcessor, ReturnValueProcessor, camel_to_snake, DatagrokFanoutType
from .func_call import FuncCall, Type, FuncCallParam, FuncCallStatus
from .settings import Settings
from .logger import get_logger
from .redirect import WebSocketRedirect, AmqpRedirect
from .amqp_publisher import AmqpFanoutPublisher


@worker_shutdown.connect
def on_worker_shutdown(**kwargs):
    if AmqpFanoutPublisher._initialized:
        try:
            AmqpFanoutPublisher.get_instance().close()
        except Exception:
            pass   

class DatagrokTask(celery.Task):
    # Disable Celery argument checking
    typing = False

    def __init__(self):
        self.settings = Settings.get_instance()
        self.logger = get_logger()

        if AmqpFanoutPublisher._initialized:
            self.amqp_publisher = AmqpFanoutPublisher.get_instance()
        else:
            self.amqp_publisher = AmqpFanoutPublisher(self.settings.broker_url, self.settings.calls_fanout)
    
        self.pipe = None

    # Overrides celery.Task.before_start
    def before_start(self, task_id, args, kwargs):
        """This function can't throw. It is special hook called before the tasks lifecycle starts, so exceptions here won't be catched 
        and task will never start.
        """
        
        self.call: FuncCall = FuncCall(args[0])
        self.accepted_send = False 
        # Try to establish AMQP connection to notify that call was accepted
        try:
            self.amqp_publisher.publish({}, self.call.id, DatagrokFanoutType.ACCEPTED)
            self.accepted_send = True
        except Exception as e:
            self.logger.warning("Initial Accepted message send failed: %s", str(e), extra={"task_id": self.call.id}, exc_info=True)   

        if self.call.requires_pipe:
            try:
                self._ensure_pipe_conn(timeout=30)     
            except Exception as e:
                self.pipe = None
                self.logger.warning("Initial grok_pipe connection failed with exception: %s", str(e), exc_info=True, extra={"task_id": self.call.id})


    # Overrides celery.Task.after_return
    def after_return(self, status, retval, task_id, args, kwargs, einfo):
        self.call.status = FuncCallStatus.COMPLETED if status == 'SUCCESS' else FuncCallStatus.ERROR

        if einfo is not None:
            self.call.error_message = str(einfo.exception)
            self.call.error_stack_trace = "".join(einfo.traceback)
        elif isinstance(retval, Exception):
            self.call.error_message = str(retval)
            self.call.error_stack_trace = ''.join(traceback.format_exception(type(retval), retval, retval.__traceback__))

        for i in range(3):
            if i > 0:
                self.logger.info("Retrying (%d) CALL send", i + 1, extra={"task_id": self.call.id})
            try:
                if self.call.requires_pipe:
                    self._pipe_send(f"CALL {json.dumps(self.call.to_json())}") 
                else:
                    self.amqp_publisher.publish(self.call.to_json(), self.call.id, DatagrokFanoutType.CALL)
                break    
            except Exception as e:
                self.logger.error("CALL send failed (attempt %d): %s", i + 1, str(e), exc_info=True, extra={"task_id": self.call.id})
                time_lib.sleep(2 ** i) 

        if self.pipe is not None:
            try:
                self.pipe.close()
            except Exception as e:
                self.logger.error("Couldn't close WS connection: %s", str(e), extra={"task_id": self.call.id}, exc_info=True)
            finally:
                self.pipe = None        

    # Overrides celery.Task.update_state
    def update_state(self, task_id=None, state=None, meta=None):
        try:
            percent = meta.get("percent") if meta else None
            description = meta.get("description") if meta else None
            if percent is None:
                self.logger.warning("Incorrect usage of update_state, meta['percent'] key should be present.")
                return
            message = {'percent': percent, 'description': description}
            if self.call.requires_pipe:
                self._pipe_send(f"PROGRESS {json.dumps(message)}")
            else:
                self.amqp_publisher.publish(message, self.call.id, DatagrokFanoutType.PROGRESS)
        except Exception as e:
            # Ignoring update status errors
            self.logger.warning("Update status failed: %s", str(e), extra={"task_id": self.call.id}, exc_info=True)           

    # Overrides celery.Task.__call__
    def __call__(self, *args, **kwargs):
        if len(self.call.output_params) > 1:
            raise ValueError("Only one return parameter is allowed, but more are present in the Datagrok function declaration.") 
         
        if not self.accepted_send:
            try:
                self.amqp_publisher.publish({}, self.call.id, DatagrokFanoutType.ACCEPTED)
                self.accepted_send = True
            except Exception as e:
                self.logger.error("Accepted message send failed: %s", str(e), extra={"task_id": self.call.id}, exc_info=True)
                raise 
        if self.call.requires_pipe:
            self._ensure_pipe_conn(timeout=30)

        self._process_input_params()

        if self.call.requires_pipe:
            with WebSocketRedirect(self.pipe):
                return self._run_and_process()
        else:
            with AmqpRedirect(self.call.id):
                return self._run_and_process()
   
    def _run_and_process(self):
        value = self._run()
        self._process_return_value(value)
        return None

    def _run(self) -> Any:
        args = []
        kwargs = {}
        param_map = {camel_to_snake(x.name): x for x in self.call.input_params}
        fn = getattr(self, '__wrapped__', None)
        if fn is not None:
            sig = inspect.signature(fn)
            parameters = [x for x in sig.parameters.values()]
            prepend_self = (
                parameters and
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
                        self.logger.warning("Celery task contains *args, skipping it", extra={"task_id": self.call.id})
                        continue
                    elif param.kind == inspect.Parameter.VAR_KEYWORD:
                        # add session token that can be used in datagrok_api.DatagrokClient
                        param_map["USER_API_KEY"] = self.call.user_api_key
                        kwargs.update(param_map)    
                        break    
            
        try:
            sig.bind(*args, **kwargs)
        except TypeError as e:
            self.logger.error("Couldn't match Datagrok func parameters to the task: %s", str(e), extra={"task_id": self.call.id}, exc_info=True)
            raise

        return self.run(*args, **kwargs)
        
    def _process_return_value(self, value):
        outputs = self.call.output_params
        if not outputs:
            return

        param = outputs[0]

        if value is None:
            param.value = None
            self.logger.warning("Task returned None", extra={"task_id": self.call.id})
            return

        setter = ReturnValueProcessor.type_map.get(param.property_type)

        if setter is None:
            raise ValueError(f"Unsupported return type: {param.property_type}")
        
        if param.property_type == Type.DATA_FRAME:
            setter(param, value, self.call.use_parquet_transfer)
        else:
            setter(param, value)

        if param.is_streamable:
            self._send_param_grok_pipe(param)

    def _send_param_grok_pipe(self, param: FuncCallParam):
        self.logger.info("Sending param %s", param.name, extra={"task_id": self.call.id})
        if not isinstance(param.value, (bytes, bytearray)):
                raise TypeError(f"Expected bytes for param.value, got {type(param.value)}")
        
        param_type = 'blob'
        param_id = param.name
        data = param.value
        size = len(data)
        if param.property_type == Type.DATA_FRAME:
            param_type = "parquet" if self.call.use_parquet_transfer else "csv"
            param_id = uuid.uuid1()
        tags_json = json.dumps({".id": param_id, ".type": param_type})

        start = time_lib.monotonic()

        self._pipe_send(f"SENDING DATAFRAME {size} {tags_json}")

        batch_size = self.call.binary_batch_size
        for i in range(0, size, batch_size):
            chunk = data[i:i+batch_size]
            self.logger.debug("Sending binary batch %s of param %s", i, param.name, extra={"task_id": self.call.id})
            self._pipe_send(chunk)
            
            try:
                response = self.pipe.recv()
            except websocket.WebSocketTimeoutException as e:
                self.logger.error("Timeout waiting for message from grok_pipe: %s", str(e), exc_info=True, extra={"task_id": self.call.id})
                raise
            if response == "ERROR":
                raise Exception(f"Error sending parameter {param.name}")
            elif response != "PART_OK":
                raise ValueError("Unexpected message received from grok_pipe.")
            
        self.logger.info("Finished sending param %s (%d bytes in %d chunks), elapsed: %d ms", param.name, size, (size + batch_size - 1) // batch_size, (time_lib.monotonic() - start) * 1000, extra={"task_id": self.call.id})
        param.value = {"id": param_id} if param.property_type == Type.DATA_FRAME else param_id    

    def _process_input_params(self):
        for param in self.call.input_params:
            if param.is_streamable:
                param.value = self._get_param_grok_pipe(param)
            if param.value is not None:
                setter = InputValueProcessor.type_map.get(param.property_type)
                if setter is not None:
                    if param.property_type == Type.DATA_FRAME:
                        setter(param, self.call.use_parquet_transfer)
                    else:
                        setter(param)
                  
    def _get_param_grok_pipe(self, param: FuncCallParam):
        self.logger.info("Receiving param %s", param.name, extra={"task_id": self.call.id})

        start = time_lib.monotonic()

        self._pipe_send(f"PARAM {param.name}")
        
        chunks = []
        total_bytes_received = 0
        expected_size = None
        
        while True:
            if (time_lib.monotonic() - start) > self.settings.param_timeout_minutes * 60:
                self.logger.error("Timeout receiving param %s", param.name, extra={"task_id": self.call.id})
                raise TimeoutError(f"Timeout receiving param {param.name}")

            try:
                message = self.pipe.recv()
            except websocket.WebSocketTimeoutException as e:
                self.logger.error("Timeout waiting for message from grok_pipe: %s", str(e), exc_info=True, extra={"task_id": self.call.id})
                raise
                
            if isinstance(message, str):
                if message.startswith(f"PARAM_SENT {param.name}"):
                    self.logger.debug("Received all chunks of %s parameter", param.name, extra={"task_id": self.call.id})
                    break
                elif message.startswith("SENDING"):
                    self.logger.debug("Server is sending parameter data: %s", message, extra={"task_id": self.call.id})
                    expected_size = self._get_expected_size(message)
                else:
                    self.logger.debug(f"Unexpected message: %s", message[:30])
            else:
                if expected_size is None:
                    self.logger.error("Expected size was not set before binary data began for param %s.", param.name, extra={"task_id": self.call.id})
                    raise ValueError("Protocol error: expected size not defined.")
                
                chunk_size = len(message)
                self.logger.debug("Received binary chunk of size %s", chunk_size, extra={"task_id": self.call.id})
                chunks.append(message)
                total_bytes_received += chunk_size
                self.logger.debug("Total bytes received: %s/%s", total_bytes_received, expected_size, extra={"task_id": self.call.id})
                self._pipe_send("PART_OK")
                if total_bytes_received > expected_size:
                    self.logger.error("Received size of binary data is larger than expected for param %s", param.name, extra={"task_id": self.call.id})
                    raise ValueError("Received size of binary data is larger than expected.")      
                if total_bytes_received == expected_size:
                    self.logger.info("Received all bytes for param %s, elapsed: %d ms", param.name, (time_lib.monotonic() - start) * 1000, extra={"task_id": self.call.id})
                    break
        
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
                expected_size = int(parts[2])
                self.logger.debug("Expected parameter size: %s bytes", expected_size)
        except (ValueError, IndexError) as e:
            self.logger.error("Could not parse size from SENDING: %s", str(e), exc_info=True, extra={"task_id": self.call.id})
            raise

    def _ensure_pipe_conn(self, max_retries=5, retry_interval=3, timeout=15):
        if self.pipe is not None and self.pipe.connected:
            return
        attempt = 0
        while attempt < max_retries:
            try:
                self.pipe = websocket.create_connection(f"{self.settings.pipe_url}/{self.call.id}", header=[f"x-member-name: celery-{self.settings.celery_name}", f"authorization: {self.settings.pipe_key}"] or [], timeout=timeout)
                self.pipe.settimeout(self.settings.ws_message_timeout_seconds)
            except Exception as e:
                self.logger.warning("WS connection failed (attempt %s/%s): %s", attempt + 1, max_retries, str(e), exc_info=True, extra={"task_id": self.call.id})
                attempt += 1
                time_lib.sleep(retry_interval)
        raise ConnectionError(f"Failed to connect to grok_pipe at {self.settings.pipe_url} after {max_retries} retries") 

    def _pipe_send(self, message):
        self._ensure_pipe_conn()
        try:
            self.pipe.send(message)
        except websocket.WebSocketConnectionClosedException:
            self.logger.warning("Grok_Pipe is disconnected, trying to reconnect...", self.call.id, extra={"task_id": self.call.id})
            self._ensure_pipe_conn()
            self.pipe.send(message) 
            