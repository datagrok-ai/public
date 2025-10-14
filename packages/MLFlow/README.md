# MLflow

`MLflow` is a [package](https://datagrok.ai/help/develop/develop#packages) for the [Datagrok](https://datagrok.ai) platform.

MLflow integrates [MLflow](http://mlflow.org) models into the Datagrok platform. It lets users fetch, manage, and run MLflow-registered models within Datagrok’s scalable infrastructure using Docker containers.

Key features:

* Model fetching: Easily access MLflow models using the MLflow connector
* Containerized execution: Deploy models in isolated Docker containers for secure and reproducible execution
* Inference API: Perform real-time predictions and visualize results within Datagrok’s interactive environment
* Version control: Access specific or latest versions of models with MLflow’s model versioning
* Environment consistency: Maintain consistent environments using MLflow’s conda or Docker specifications
* Seamless integration: Combine Datagrok’s analytics and visualization with ML model inference

Use cases:

* Perform automated model inference on large datasets within Datagrok
* Collaborate on data science workflows integrating data wrangling, visualization, and machine learning