# Datagrok Python Celery Task Library

This package provides a Celery task wrapper for integrating Python functions into the [Datagrok](https://datagrok.ai) platform. It enables logging, progress tracking, and seamless communication between Celery workers and Datagrok’s infrastructure.

Refer to the [Datagrok Help](https://datagrok.ai/help/datagrok) for more information about the platform.

---

## Installation

Install the package using [pip](https://pypi.org/project/pip/):

```bash
pip install datagrok-celery-task
```

## Usage

To define tasks compatible with Datagrok, use the DatagrokTask base class and configure your Celery app with the provided Settings class.

```python
from celery import Celery
from datagrok_celery_task import DatagrokTask, Settings
import logging

# Always create a Settings object. Provide properties manually only if not launched by Datagrok.
settings = Settings(log_level=logging.DEBUG)

# Create a Celery app
app = Celery(settings.celery_name, broker=settings.broker_url)

# Define a simple Datagrok task
@app.task(base=DatagrokTask)
def echo(c, **kwargs):
    print("Received USER_API_KEY:", kwargs.get("USER_API_KEY", "empty"))
    print("Received data:", c)

# Define a task that reports progress
@app.task(bind=True, base=DatagrokTask)
def progress_task(self: DatagrokTask, a):
    self.update_state(meta={"percent": 10, "description": "Starting"})
    print(a)
    return "Task completed"
```

### Notes

* When Datagrok manages the Celery worker, environment variables will auto-populate Settings.
* Prins will appear in Datagrok log messages of the FuncCall. Please note, that it will work only when Celery `prefork worker pool` is used.
* Bind task `bind=True` and use `self.update_state` to manually update platform's progress bar.
* Tasks can return results or raise Exception.

## License

See [License.md](https://github.com/datagrok-ai/public/blob/master/LICENSE.md).
