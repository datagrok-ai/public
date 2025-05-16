# Datagrok python client library

This library can be used for integration with datagrok. It is a Celery task wrapper to wrap python funcitons

Refer to [Help](http://datagrok.ai/help) for more information.

## Installation

To install package, use [pip](https://pypi.org/project/pip/).

```shell
pip install datagrok-celery-task
```

## Usage

To use API client, import DatagrokClient to your project:

```python
from celery import Celery
from datagrok_celery_task import DatagrokTask, Settings
import logging

# Init it with correct properties if celery is started not by the Datagrok infrastructure.
# If Datagrok starts it then Settings properties will be populated from env variables.
settings = Settings(log_level=logging.DEBUG)
app = Celery(settings.celery_name, broker=settings.broker_url)

# Use DatagrokTask as base for your tasks
@app.task(base=DatagrokTask)
def test(c, **kwargs):
    # All print statements will be sent to Datagrok platform as log messages
    print("Got kwargs: " + kwargs.get("USER_API_KEY", "empty"))
    print("Got data:" + str(c))

# Add bind flag for your function to be able to control progress indication
@app.task(bind=True, base=DatagrokTask)
def debug_task(self: DatagrokTask, a):
    # All print statements will be sent to Datagrok platform as log messages
    self.update_state(meta={"percent": 10, "description": "Running something"})
    print(a)
    return "Hello World"
```

## License

See [License.md](https://github.com/datagrok-ai/public/blob/master/LICENSE.md).
