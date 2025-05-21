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