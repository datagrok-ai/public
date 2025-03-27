import os
import logging
import asyncio
from quart import Quart, request
from quart import websocket as quart_websocket
from hypercorn.asyncio import serve
from hypercorn.config import Config

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Quart app (ASGI-compatible Flask alternative)
app = Quart('apitests')

@app.route('/square', methods=['GET'])
async def square():
    num = int(request.args.get("number", 5))
    logger.info(f"Calculating square of {num}")
    return {"result": num * num}, 200

@app.route('/health', methods=['GET'])
async def health():
    num = 2 * 2
    if num == 4:
        logger.info("Health check: Healthy")
        return {"result": "Healthy"}, 200
    else:
        logger.error("Health check: Unhealthy")
        return {"result": "Unhealthy"}, 500

# WebSocket handler
@app.websocket('/ws')
async def websocket_handler():
    logger.info("WebSocket connection established")
    while True:
        message = await quart_websocket.receive()
        logger.info(f"Received message: {message}")
        await quart_websocket.send(message)

# Run the ASGI server
async def run_server():
    config = Config()
    config.bind = ["0.0.0.0:5353"]  # Same port for HTTP and WebSocket
    await serve(app, config)

if __name__ == '__main__':
    logger.info("Starting ASGI server on port 5353")
    asyncio.run(run_server())