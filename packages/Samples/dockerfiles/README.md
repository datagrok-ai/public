# Docker for custom model

This folder contains a Docker configuration for running a [custom model](https://datagrok.ai/help/learn/custom-machine-learning-models) as a Celery task in Datagrok:

```
dockerfiles/  
  ├── Dockerfile        # Container setup and dependencies  
  ├── app.py            # Celery entrypoint  
  ├── container.json    # Container metadata  
  ├── environment.yaml  # Conda environment  
  ├── predict.py        # Model inference  
  └── train.py          # Model training  
```

For details, see the [Celery Python functions documentation](https://datagrok.ai/help/develop/how-to/packages/python-functions).