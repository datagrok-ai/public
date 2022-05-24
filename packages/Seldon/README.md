# Seldon

Seldon is a [package](https://datagrok.ai/help/develop/develop#packages) for the
[Datagrok](https://datagrok.ai) platform. It allows connecting to a
[Seldon Cluster](https://www.seldon.io/) with deployed machine learning models and easily apply them within Datagrok.

## Overview

Seldon is a popular ML Ops platform for managing and using containerized machine learning models. Many of such models
are of great use directly inside Datagrok through its built-in
[custom machine learning models](https://datagrok.ai/help/learn/custom-machine-learning-models)
infrastructure. The package `Seldon` uses this infrastructure to add support for running arbitrary Seldon models simply
by pointing to the cluster of interest. Models with their inputs and outputs need not be known in advance in order to be
used. It is also possible to enumerate and connect to Seldon models programmatically in your Datagrok extensions and
applications in TypeScript/JavaScript.

Since version 1.3 Seldon supports [metadata](https://github.com/SeldonIO/ml-prediction-schema), or prediction schema,
which allows understanding, for any given model, its inputs, outputs, their order, types, inner structure (one-hot,
probabilistic) and other traits. `Seldon` package supports any Seldon-deployed model which provides for this metadata.

Models are deployed to a Seldon Cluster into **namespaces** and are called **deployments**. In Seldon there is a
new [Seldon Deployment SDK](https://github.com/SeldonIO/seldon-deploy-sdk)
which provides for all the relevant enumeration of these namespaces and deployments (models).

Currently, the package support identification with [Keycloak](https://www.keycloak.org/).

## Usage

After deploying the package `Seldon`, set up the followng [credentials]() for it
(the values are provided as examples):

* seldonUser
* seldonPassword
* seldonHost: `https://models.YOUR-SERVER.com/seldon-deploy/api/v1alpha1`
* seldonOIDCServer: `https://isaac-keycloak.YOUR-SERVER.com/auth/realms/orgName`
* seldonClientID: `seldon-api`
* seldonNamespace: `test-namespace`

## Roadmap

### Alpha version

Is built for pre-v1.3 deployments which don't have model metadata provided. The demo mocks the metadata for a
popular `sklearn` predictive model `iris`.

Open [`iris.csv`]() in Datagrok, activate the [top menu]() and find the item `Seldon | Apply`. The list of deployments (
models) will be loaded for a pre-specified namespace. Select the one corresponding to `iris` and hit `Ok`. The `iris`
dataframe will be shortly filled in with predictions coming from the Seldon model.

### Beta version

* Supports enumeration of the namespaces
* Supports arbitrary metadata parsing for Seldon versions >= 1.3

### Release 1.0

* Supports all kinds of data types supported by Seldon metadata
