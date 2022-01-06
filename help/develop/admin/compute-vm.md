<!-- TITLE: Compute VM -->
<!-- SUBTITLE: -->

# Compute VM

Compute VM is a set of [docker containers](https://www.docker.com/) that includes all Datagrok platform computational
services. Server is easy scalable using [Amazon ECS](https://aws.amazon.com/ecs/)
or [Kubernetes](https://kubernetes.io/) services.

## Services

![Compute VM](../../uploads/features/compute-vm.png "Compute VM")

| Service                                                             | Docker Image | Version  | Port  |
|---------------------------------------------------------------------|--------------|----------|-------|
| [GrokCompute](grok-compute.md)                                      | [latest](https://hub.docker.com/r/datagrok/grok_compute) | 1.3.0    | 5005  |
| [GrokHelper](architecture-details.md#grok-helper)                   | - | 1.0.0    | 5005  |
| [Jupyter Notebooks](https://jupyter.org)                            | [latest](https://hub.docker.com/r/datagrok/jupyter_notebook) | 6.0.3    | 8889  |
| [Jupyter Kernel Gateway](https://github.com/jupyter/kernel_gateway) | [latest](https://hub.docker.com/r/datagrok/jupyter_kernel_gateway) | 2.5.1    | 8888  |
| [H2O](https://www.h2o.ai/products/h2o/)                             | [latest](https://hub.docker.com/r/datagrok/h2o) | 3.26.0.5 | 54321 |

## Programming languages

| Service                                      | Version  |
|----------------------------------------------|----------|
| [Python](https://www.python.org)             | 3.6.15   |
| [R](https://www.r-project.org)               | 4.1.2    |
| [Julia](https://julialang.org)               | 1.6.2    |
| [NodeJS](https://nodejs.org)                 | 16.8.0   |
| [Octave](https://octave.sourceforge.io/)     | 5.2.0    |

## Installation

* [Compute VM regular deploy](deploy-regular.md#setup-compute-virtual-machine)
* [Compute VM Amazon EC2 deploy](deploy-amazon-ec2.md#setup-compute-virtual-machine)
* [Compute VM Amazon ECS deploy](deploy-amazon-ecs.md#setup-compute-virtual-machine)

## Scalability

Compute VM is easy scalable in horizontal direction. For example: in case of setup in
[AWS ECS](https://aws.amazon.com/ecs/)
cluster [Application Load Balancer](https://aws.amazon.com/elasticloadbalancing/)
can be used.

At the moment some limitations are exists:

* "Substructure search" and "Descriptors" API in [GrokCompute](grok-compute.md) should not be used in cache mode.

See also:

* [Compute Virtual Machine Architecture](architecture-details.md#compute-virtual-machine)
* [GrokCompute](grok-compute.md)
