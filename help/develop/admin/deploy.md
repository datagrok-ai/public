<!-- TITLE: Deployment -->
<!-- SUBTITLE: -->

# Deployment

This document contains instructions to deploy Datagrok.

## Prerequisites

1. Provide required [resources for every component](infrastructure.md#resources)
2. Both [Compute](infrastructure.md#compute-components) and [Datagrok](infrastructure.md#datagrok-components) engines
   should be accessible by users. The easiest way is to create DNS endpoints pointing to load balancers in front of the
   services:
   `datagrok.example` and `cvm.example`.
3. Configure [database](infrastructure.md#database) and [storage](infrastructure.md#storage), which should be available
   for [Datlas](infrastructure.md#datlas). If you use local setup with local file storage and local database, skip this
   step.

## 