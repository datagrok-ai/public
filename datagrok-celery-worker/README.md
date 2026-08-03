# datagrok-celery-worker

A generic Node.js, Celery-compatible worker that executes a published Datagrok package's
JS functions. It is the JS mirror of the python library
[`datagrok-celery-task`](../datagrok-celery-task): the Datagrok server (datlas) queues a
`DockerFunc` call into a RabbitMQ task queue, the worker picks it up, streams dataframe /
blob parameters over [grok_pipe](../../grok_pipe), runs the function from the package's
published webpack bundle (loaded via the `datagrok-api` Node runtime), and publishes the
result back.

## How it works

1. **Task pickup** — consumes the queue named by `TASK_QUEUE_NAME` (the worker owns the
   `assertQueue {durable: true}` declaration; datlas publishes into it passively). The
   message body is `{"args": [funcCallJson], "kwargs": {}}` (celery protocol-2 array
   bodies are accepted too). The message is acked immediately and an `accepted` message
   is published to the `calls_fanout` exchange — datlas requires it within its
   task-pickup timeout (60 s by default).
2. **Package loading** — on the first task the worker calls
   `startDatagrok({apiUrl, apiToken, detached: true})` and `loadPackage(name)` from
   `datagrok-api/datagrok`, then resolves the function implementation directly from the
   bundle's module exports (never through `grok.functions.call`, which would recurse to
   the server entity).
3. **Param streaming** — when the call has dataframe/blob/file params, the worker joins
   the grok_pipe room `ws://$DATAGROK_PIPE_HOST:$DATAGROK_PIPE_PORT/<callId>` (headers
   `x-member-name: celery-$DATAGROK_CELERY_NAME`, `authorization: $DATAGROK_PIPE_KEY`).
   Inputs: send `PARAM <name>`, receive `SENDING DATAFRAME <size> <tagsJson>` + binary
   chunks (each acked with `PART OK`), terminated by `PARAM_SENT <name>`. Outputs: the
   same `SENDING`/chunks/`PART OK` exchange in the other direction. Dataframes travel as
   native d42 binary (`DG.DataFrame.fromByteArray` / `df.toByteArray()`, pipe tag
   `.type: 'dataframe'`); CSV input (`.type: 'csv'`) is still accepted for compatibility
   with older datlas versions that send CSV to js-lang funcs.
4. **Result** — a `CALL <funcCallJson>` text frame over the pipe when one is open,
   otherwise a `call` message on `calls_fanout`. Progress reports go out as
   `PROGRESS {json}` frames / `progress` fanout messages; package code can call
   `globalThis.DG_TASK_PROGRESS(percent, description)`.
5. **Cancellation** — datlas publishes `{method: 'revoke', arguments: {task_id}}` to the
   `celery.pidbox` exchange; the worker keeps a per-hostname pidbox queue
   (`celery@$CELERY_HOSTNAME.celery.pidbox`). A revoke immediately publishes a
   `Canceled` result and sets a cooperative cancel flag — a running function is not
   killed, but its next progress call throws.

## Environment variables

| Variable | Required | Default | Description |
|----------|----------|---------|-------------|
| `TASK_QUEUE_NAME` | yes | | RabbitMQ task queue to consume (the DockerFunc's `task_queue`) |
| `CELERY_HOSTNAME` | yes | | Worker hostname; forms the pidbox queue name `celery@<hostname>.celery.pidbox` |
| `DATAGROK_PACKAGE_NAME` | yes | | Published Datagrok package whose JS functions this worker executes |
| `DATAGROK_PACKAGE_VERSION` | no | | Logged only; the worker always loads the server's current published version |
| `DATAGROK_CELERY_NAME` | no | `datagrok-celery` | grok_pipe member name suffix (`celery-<name>`) |
| `DATAGROK_AMPQ_HOST` | no | `localhost` | RabbitMQ host (note the existing 'AMPQ' spelling, kept from the python lib) |
| `DATAGROK_AMPQ_PORT` | no | `5672` | RabbitMQ port |
| `DATAGROK_AMPQ_USER` | no | `guest` | RabbitMQ user |
| `DATAGROK_AMPQ_PASSWORD` | no | `guest` | RabbitMQ password |
| `DATAGROK_AMPQ_TLS` | no | `false` | `true` → `amqps://` |
| `DATAGROK_PIPE_HOST` | no | `localhost` | grok_pipe host |
| `DATAGROK_PIPE_PORT` | no | `3000` | grok_pipe port |
| `DATAGROK_PIPE_KEY` | no | empty | grok_pipe `authorization` header value |
| `DATAGROK_API_URL` | no | | Fallback Datagrok API url; normally each call carries `aux.DATAGROK_API_URL` |
| `DATAGROK_PARAM_TIMEOUT` | no | `5` | Total budget for receiving one param, minutes |
| `DATAGROK_WS_MESSAGE_TIMEOUT` | no | `30` | Per-message pipe wait timeout, seconds |
| `HEALTH_PORT` | no | `8000` | Port of the `{"status":"ok"}` health endpoint |

## Local run

Start RabbitMQ and grok_pipe from the datlas compose stack:

```bash
cd core/server/datlas/resources/compose
docker compose -f deps.docker-compose.yaml -p datagrok-dev --profile scripting up -d rabbitmq grok_pipe
```

Then build and run the worker against them (with datlas running on 8082):

```bash
npm install && npm run build
TASK_QUEUE_NAME=my_pkg_queue \
CELERY_HOSTNAME=local-worker \
DATAGROK_PACKAGE_NAME=MyPackage \
DATAGROK_API_URL=http://localhost:8082/api \
npm start
```

`datagrok-api` is a peer dependency — install a version that exposes the
`datagrok-api/datagrok` Node entrypoint (>= 1.27) next to this package.

## Limitations

- **No parquet dataframe transfer** — dataframes travel as native d42 binary in both
  directions (preserving column types and tags); CSV input is still accepted from an
  older datlas. A call with `options.isParquet = true` fails fast.
- **Single in-flight task** — tasks are acked immediately and drained from an in-memory
  FIFO one at a time.
- **Cooperative cancellation** — a revoke publishes the `Canceled` result right away,
  but the running JS function is only interrupted at its next
  `DG_TASK_PROGRESS` call.
- **One return parameter max** (same rule as the python lib).

## Security

Never log `aux` or `USER_API_KEY` — the call's `aux` carries the user's session token.
The worker only echoes `aux` back inside the result message (as the protocol requires)
and never prints it.
