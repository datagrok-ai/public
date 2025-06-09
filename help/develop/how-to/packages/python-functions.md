# Creating Python Docker Apps

## Overview

This guide outlines how to build and deploy your own Python-based Docker applications that integrate seamlessly with the Datagrok platform. Your Python functions will be executed as [Celery](https://docs.celeryq.dev/en/stable/index.html#) tasks within isolated Docker containers. The platform takes care of container orchestration, task queuing, and result handling.

## Key Benefits

* **Write only the logic**: You focus on writing plain Python functions.
* **Zero infrastructure setup**: The platform handles Celery, RabbitMQ, Docker, and scaling.
* **Full control**: Optional configuration for resource limits and dependencies.

---

## Folder Structure

Create a `python/` directory inside your plugin root. Inside it:

* One or more folders, each representing a separate application
* Each folder can contain:

    * One or more `.py` files with functions.
    * An optional `requirements.in` or `environment.yaml` file for Python dependencies
    * An optional `container.json` file for resource configuration

**Example layout:**

```
plugin-root/
└── python/
    └── my_app/
        ├── logic.py
        ├── requirements.in
        └── container.json
```

---

## Writing Your Python Code

You can define any number of functions inside your Python files. Functions must include metadata as comments:

```python
# file: logic.py

#name: add
#tags: task
#input: int x
#input: int y
#output: int z
def add(x, y):
    return x + y
```

This metadata allows the platform to expose the function in the UI and JS API.

---

## Optional Files

### `requirements.in`

List of pip dependencies for your app:

```text
numpy
scikit-learn
```


### `environment.yaml`

Valid [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) if you want to use miniconda as your package manager.

```text
name: my-simple-app
channels:
  - defaults
dependencies:
  - python=3.11
  - pip
```

### `container.json`

Optional resource [configuration](docker-containers.md#2-container-configuration):

```json
{
  "cpu": 1,
  "memory": 1024
}
```

By default, [Celery prefork pool](https://docs.celeryq.dev/en/stable/internals/reference/celery.concurrency.prefork.html#module-celery.concurrency.prefork) is used and the number of process workers will be `cpu * 2`.

> Note: The default and minimum value of `cpu` parameter for Celery based containers is 1

---

## Platform Responsibilities

When you deploy a plugin:

1. **Function Parsing**: The platform parses your annotated Python functions.
2. **Task Wrapping**: A Celery-main file is generated automatically.
3. **Docker Build**: A Dockerfile is created using boilerplate.
4. **Container Deployment**: Containers are started using `grok_spawner`.
5. **Celery Worker Setup**: Each container runs a Celery worker bound to your task queue.

---

## Executing Tasks

Once your app is deployed, functions can be called via the JS API:

```javascript
await grok.functions.call('Plugin:add', { x: 1, y: 2 });
```

The platform will:

* Publish the task to RabbitMQ
* Route it to the correct worker container
* Execute your function
* Collect the results and return it back to you

---

## Best Practices

* Keep your apps stateless to ensure parallel execution safety
* Use metadata comments to describe inputs/outputs clearly
* Declare dependencies explicitly in `requirements.in`
* Test locally before deploying

---

## Troubleshooting

* **Task not showing up?** Check function annotations for correct formatting
* **Dependency error?** Validate your `requirements.in` by running `pip install -r requirements.in` locally
* **Performance issues?** Tune `cpu` and `memory` in `container.json`

---

## Conclusion

This approach allows you to define lightweight, scalable, and easily deployable Docker-based Python apps with minimal effort. By leveraging Celery and RabbitMQ, the platform ensures
high performance, reliability, and scalability out of the box.


See also:
- [Packages Docker containers](docker-containers.md)
- [Packages](../../develop.md#packages)
- [Connecting to database inside package Docker container](../db/access-data.md)
