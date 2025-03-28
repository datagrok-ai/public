import pika
import json
import os
import base64
import pandas as pd
from io import BytesIO
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['figure.figsize'] = (8.0, 8.0)
import websockets
import asyncio
import struct
from typing import Union, Dict, List, Any, Optional
from concurrent.futures import ThreadPoolExecutor
import threading
import logging
import sys


# Configure logging to ensure it flushes to Docker logs
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.StreamHandler(stream=sys.stderr)  # Output to stderr for Docker logs
    ]
)
logger = logging.getLogger("python_docker")
logger.setLevel(logging.INFO)  # Set to INFO, DEBUG, etc as needed

# Force stdout/stderr to be unbuffered (helps with Docker logs)
sys.stdout.reconfigure(line_buffering=True) 
sys.stderr.reconfigure(line_buffering=True)


def log_info(message):
    """Log info message and force flush"""
    logger.info(message)
    sys.stderr.flush()  # Force flush to Docker logs


def log_error(message):
    """Log error message and force flush"""
    logger.error(message)
    sys.stderr.flush()  # Force flush to Docker logs


def log_debug(message):
    """Log debug message and force flush"""
    logger.debug(message)
    sys.stderr.flush()  # Force flush to Docker logs


def convert_graphics_to_b64(image):
    plt.tight_layout()
    buffered = BytesIO()
    image.savefig(buffered, format="PNG", bbox_inches='tight', pad_inches=0.1)
    buffered.seek(0)
    return base64.b64encode(buffered.getvalue()).decode('utf-8')

def convert_pandas_to_csv(df):
    return df.to_csv().encode('utf-8')


def execute_code(code, exec_locals, res_var = None, res_type = None):
    response = {"status": "ok", "output": None, "error": None}
    try:
        exec_globals = {}
        exec(code, exec_globals, exec_locals)
        returned = exec_locals.get(res_var, f"No '{res_var}' variable returned")
        if res_type == 'graphics':
            returned = convert_graphics_to_b64(returned)
        elif res_type == 'dataframe':
            returned = convert_pandas_to_csv(returned)
        response["output"] = returned
    except Exception as e:
        response["status"] = "error"
        response["error"] = str(e)
    return response


async def get_script_param(websocket, param_name: str, timeout=30, debug_logger=log_debug):
    """
    Retrieves a parameter from the websocket pipe following the protocol in scripting.dart.
    
    Args:
        websocket: The websocket connection
        param_name: Name of the parameter to retrieve
        timeout: Maximum seconds to wait for completion
        debug_logger: Optional logger for debugging
        
    Returns:
        The parameter data as bytes or deserialized object
    """
    if debug_logger:
        debug_logger(f"Requesting parameter: {param_name}")
    
    # Request the parameter
    await websocket.send(f"PARAM {param_name}")
    
    # Collect chunks of data
    chunks = []
    total_bytes_received = 0
    expected_size = None
    start_time = asyncio.get_event_loop().time()
    
    try:
        while True:
            # Check for timeout
            if asyncio.get_event_loop().time() - start_time > timeout:
                debug_logger(f"Timeout waiting for parameter {param_name} after {timeout} seconds")
                if chunks:
                    debug_logger(f"Got {len(chunks)} chunks before timeout, will use them anyway")
                    break
                raise TimeoutError(f"Timed out waiting for parameter {param_name}")
            
            # Receive message with shorter timeout
            try:
                message = await asyncio.wait_for(websocket.recv(), timeout=5.0)
            except asyncio.TimeoutError:
                debug_logger("No message received for 5 seconds, continuing to wait...")
                continue
                
            # Process the message
            if isinstance(message, str):
                debug_logger(f"Received string from websocket: {message[:100]}")
                if message.startswith(f"PARAM_SENT {param_name}"):
                    debug_logger(f"Received all chunks of {param_name} parameter (PARAM_SENT)")
                    break
                elif message.startswith("SENDING"):
                    debug_logger(f"Server is sending parameter data: {message}")
                    # Try to parse the size from the message
                    try:
                        # Format: "SENDING DATAFRAME 4355 {json}"
                        parts = message.split()
                        if len(parts) >= 3:
                            expected_size = int(parts[2])
                            debug_logger(f"Expected parameter size: {expected_size} bytes")
                    except (ValueError, IndexError) as e:
                        debug_logger(f"Could not parse size from SENDING message: {e}")
                else:
                    debug_logger(f"Unexpected message: {message}")
            else:
                # Binary data chunk
                chunk_size = len(message)
                debug_logger(f"Received binary chunk of size {chunk_size}")
                chunks.append(message)
                total_bytes_received += chunk_size
                debug_logger(f"Total bytes received: {total_bytes_received}/{expected_size if expected_size else 'unknown'}")
                
                # Send acknowledgment
                await websocket.send("PART_OK")
                
                # Reset the timeout timer after receiving data
                start_time = asyncio.get_event_loop().time()
                
                # If we know the expected size and have received it all, break
                if expected_size is not None and total_bytes_received >= expected_size:
                    debug_logger(f"Received all expected bytes ({total_bytes_received}/{expected_size}), exiting loop")
                    break
        
        # Combine all chunks
        if len(chunks) == 1:
            debug_logger(f"Returning single chunk of size {len(chunks[0])}")
            return chunks[0]
        else:
            combined = b''.join(chunks)
            debug_logger(f"Returning combined chunks of total size {len(combined)}")
            return combined
    except Exception as e:
        debug_logger(f"Error in get_script_param: {str(e)}")
        raise


async def post_script_param(websocket, param_name: str, data: Union[bytes, bytearray], 
                          param_type: str = None, param_id: str = None, 
                          batch_size: int = 2048000, debug_logger=log_debug):
    """
    Sends a parameter through the websocket pipe following the protocol in scripting.dart.
    
    Args:
        websocket: The websocket connection
        param_name: Name of the parameter to send
        data: The binary data to send
        param_type: Type of parameter ('csv', 'parquet', 'blob')
        param_id: ID of the parameter value
        batch_size: Size of each chunk to send
        debug_logger: Optional logger for debugging
    """
    if not isinstance(data, (bytes, bytearray)):
        raise ValueError("Data must be bytes or bytearray")
    
    # Default to BLOB type if not specified
    if not param_type:
        param_type = "blob"
        
    # Default param_id to the param name for blob types
    if not param_id:
        param_id = param_name
    
    # Format tags according to expected format: {".id": "paramId", ".type": "type"}
    tags_json = json.dumps({".id": param_id, ".type": param_type})
    
    size = len(data)
    
    if debug_logger:
        debug_logger(f"Sending parameter {param_name} with {size} bytes, type: {param_type}, id: {param_id}")
    
    try:
        # Notify that we're starting to send data
        await websocket.send(f"SENDING DATAFRAME {size} {tags_json}")
        
        # Send data in chunks
        for i in range(0, size, batch_size):
            chunk = data[i:i+batch_size]
            if debug_logger:
                debug_logger(f"Sending chunk {i//batch_size + 1}/{(size+batch_size-1)//batch_size}, size: {len(chunk)}")
            await websocket.send(chunk)
            
            # Wait for acknowledgment
            response = await websocket.recv()
            if response == "ERROR":
                raise Exception(f"Error sending parameter {param_name}")
            elif response != "PART_OK":
                if debug_logger:
                    debug_logger(f"Unexpected response: {response}")
        
        if debug_logger:
            debug_logger(f"Sent parameter {param_name} successfully")
    except Exception as e:
        if debug_logger:
            debug_logger(f"Error in post_script_param: {str(e)}")
        raise


# Event loop for each thread
thread_local = threading.local()

def get_event_loop():
    """Get or create an event loop for the current thread."""
    if not hasattr(thread_local, "loop"):
        thread_local.loop = asyncio.new_event_loop()
        asyncio.set_event_loop(thread_local.loop)
    return thread_local.loop


async def process_request_async(uri, request):
    """Process a request asynchronously."""
    websocket = None
    try:
        code = request.get("func", {}).get('script', '')
        params = request.get('func', {}).get('params', [])
        outputParam = None
        outputType = None
        inputs = {}
        inputs_types = {}
        
        for param in params:
            if param.get('isInput', True) == False:
                outputParam = param['name']
                outputType = param['propertyType']
            else:
                inputs[param['name']] = None
                inputs_types[param['name']] = param['propertyType']
        
        for param, paramValue in request.get('parameterValues', {}).items():
            if paramValue is not None:
                inputs[param] = paramValue
        
        log_info(f"Connecting to websocket at {uri}")
        # Create websocket connection for param handling
        async with websockets.connect(uri, 
                                       additional_headers={"Authorization": 'test-key', 'x-member-name': 'python_docker'},
                                       ping_interval=10,  # Keep connection alive
                                       ping_timeout=30,   # Wait longer for pong
                                       close_timeout=5) as websocket:  # Wait for close frame
            log_info("WebSocket connected")
            
            # Handle params that need to be retrieved from the websocket
            for param, paramType in inputs_types.items():
                if paramType in ['blob', 'dataframe', 'file'] and param in inputs:
                    log_info(f"Getting parameter {param} of type {paramType}")
                    try:
                        inputs[param] = await get_script_param(websocket, param)
                        log_info(f"Successfully retrieved {param}, got data of size: {len(inputs[param]) if inputs[param] is not None else 'None'}")
                        if paramType == 'dataframe':
                            inputs[param] = pd.read_csv(BytesIO(inputs[param]))
                    except Exception as e:
                        log_error(f"Error getting parameter {param}: {str(e)}")
                        raise
            
            # Execute the code
            log_info("Executing code")
            response = execute_code(code, inputs, outputParam, outputType)

            if response['status'] == 'ok':
                request['parameterValues'][outputParam] = response['output']
                request['status'] = 'Completed'

                if outputType == 'graphics':
                    aux_map = request.get('aux', {})
                    aux_map[outputParam] = 'image/png'
                    request['aux'] = aux_map
            else:
                request['status'] = 'Error'
                request['errorMessage'] = str(response['error'])
            
            # Prepare to send the result
            result = request
            
            # Handle sending large output parameters via the param protocol
            if outputType in ['file', 'dataframe', 'blob']:
                output_data = response['output']
                # For dataframe, use CSV or PARQUET format
                param_id = outputParam
                param_type = 'blob'
                if outputType == 'dataframe':
                    use_parquet = result.get('options', {}).get('IS_PARQUET_KEY', False)
                    param_type = "parquet" if use_parquet else "csv"
                    for param in params:
                        if param['name'] == outputParam:
                            param_id = param['id']
                
                log_info(f"Sending output parameter {outputParam}")
                result['parameterValues'][outputParam] = param_id
                try:
                    await post_script_param(
                        websocket, 
                        outputParam, 
                        output_data, 
                        param_type=param_type,
                        param_id=param_id
                    )
                    log_info(f"Successfully sent output parameter {outputParam}")
                except Exception as e:
                    log_error(f"Error sending parameter: {str(e)}")
                    raise       

            # Send response over websocket
            log_info("Sending CALL message with results")
            await websocket.send('CALL ' + json.dumps(result))
            log_info('SENT CALL message')

            return result
            
    except Exception as e:
        log_error(f"Error in process_request_async: {str(e)}")
        if websocket is not None and not websocket.closed:
            try:
                await websocket.close()
                log_info("Forcibly closed WebSocket connection after error")
            except:
                pass
        request['status'] = 'Error'
        request['errorMessage'] = str(e)
        return request


def on_request(ch, method, properties, body):
    """Synchronous handler for RabbitMQ messages that dispatches async work to the event loop."""

    request = json.loads(body)
    ch.basic_ack(delivery_tag=method.delivery_tag)
    
    if request.get('status', 'Completed') == 'Completed':
        return
    
    # Acknowledge that we've accepted the task
    ch.basic_publish(
        exchange='calls_fanout',
        routing_key='',
        properties=pika.BasicProperties(type='accepted', correlation_id=properties.correlation_id),
        body=''
    )

    log_info(f"Processing request with correlation ID: {properties.correlation_id}")
    
    # Setup the WebSocket URI
    uri = f'{os.environ['DATAGROK_PIPE_HOST']}:{os.environ['DATAGROK_PIPE_PORT']}/{properties.correlation_id}'
    
    # Get the event loop for this thread
    loop = get_event_loop()
    
    # Run the async processing in the event loop
    try:
        response = loop.run_until_complete(process_request_async(uri, request))
    except Exception as e:
        log_error(f"Error in async processing: {str(e)}")
        request['status'] = 'Error'
        request['errorMessage'] = str(e)
        response = request
    
    log_info(f"Finished processing request, status: {response.get('status')}")


def main():
    log_info("===== STARTING PYTHON DOCKER SERVER =====")
    try:
        # Print environment variables (without sensitive info)
        log_info("Environment setup:")
        log_info(f"- PYTHONUNBUFFERED: {os.environ.get('PYTHONUNBUFFERED', 'not set')}")
        
        amqp_host = os.environ['DATAGROK_AMPQ_HOST']
        amqp_port = os.environ['DATAGROK_AMPQ_PORT']
        amqp_user = os.environ['DATAGROK_AMPQ_USER']
        amqp_password = "********"  # Don't log the actual password
        
        log_info(f"Connecting to RabbitMQ at {amqp_host}:{amqp_port}")
        
        credentials = pika.PlainCredentials(os.environ['DATAGROK_AMPQ_USER'], os.environ['DATAGROK_AMPQ_PASSWORD'])
        parameters = pika.ConnectionParameters(
            host=amqp_host, 
            port=int(amqp_port), 
            credentials=credentials,
            heartbeat=600,  # Increase heartbeat to prevent connection drops
            blocked_connection_timeout=300
        )
        connection = pika.BlockingConnection(parameters)
        channel = connection.channel()
        queue_name = 'python_docker'
        log_info(f"Deleting existing queue: {queue_name}")
        channel.queue_delete(queue=queue_name)
        log_info(f"Creating queue: {queue_name}")
        channel.queue_declare(queue=queue_name)
        
        log_info(f"Setting up consumer for queue: {queue_name}")
        channel.basic_consume(queue=queue_name, on_message_callback=on_request)
        
        log_info("===== READY: Waiting for code execution requests... =====")
        # Flush all logs before start_consuming
        sys.stdout.flush()
        sys.stderr.flush()
        
        channel.start_consuming()
    except KeyboardInterrupt:
        log_info("Stopping server due to keyboard interrupt...")
    except Exception as e:
        log_error(f"Error in main: {str(e)}")
        import traceback
        log_error(traceback.format_exc())


if __name__ == "__main__":
    # Make sure Python doesn't buffer stdout/stderr
    os.environ['PYTHONUNBUFFERED'] = '1'
    main()
