import os
import importlib
import logging

import yaml
from celery import Celery
from datagrok_celery_task import DatagrokTask, Settings

import sys
sys.path.insert(0, os.getcwd())

settings = Settings(log_level=logging.DEBUG)
app = Celery(settings.celery_name, broker=settings.broker_url)

with open("tasks.yaml") as f:
    tasks = yaml.safe_load(f)

if not isinstance(tasks, dict) or not isinstance(tasks['tasks'], list):
    raise RuntimeError('Incorrect format for tasks.yaml. It should contain tasks field.')

for task_def in tasks['tasks']:
    mod = importlib.import_module(task_def["module"])
    fn = getattr(mod, task_def["name"])
    app.task(name=task_def["name"], base=DatagrokTask)(fn)
