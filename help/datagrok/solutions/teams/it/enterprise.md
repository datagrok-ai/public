---
title: "Datagrok for Enterprise IT"
description: Overview of enterprise IT capabilities in Datagrok, including security, storage, governance, and monitoring.
keywords:
  - role-based access control
  - authentication and authorization
  - data governance
  - on-premises deployment
  - usage monitoring
  - audit trail
  - grok compute
---

Software is eating the world, and IT has become a critical competence in all businesses. Datagrok enables IT to become a
partner with business by solving business problems, while at the same time providing enterprise-grade security,
governance, and other IT-controlled mechanisms for storage administration, backup, and audit.

## Security

Give the right access to the right people. We offer several
[authentication](../../../../govern/access-control/access-control.md#authentication) and [authorization](../../../../govern/access-control/access-control.md#authorization) options, as well as role- and group-based
privileges for all objects available in the platform.

[Learn more](../../../../govern/access-control/access-control.md).

## Storage

You are in full control of your data. Pick whatever storage your organization prefers, or mix and match. We support
network file systems, S3, Dropbox, and Google Cloud.

## Governance

Centrally manage all of the data sources, file shares, connections, queries, and reports in one place. See which
queries, scripts, and models are behind a result with [data provenance](../../../../govern/audit/data-provenance.md),
and who uses which data with [Usage Analysis](../../../../govern/audit/usage-analysis.md).

## Deployment

We are completely platform agnostic, you can deploy Datagrok on Linux, Windows, on-premises, in the private cloud, or
use it as a SaaS. To scale scientific computations, spin out as many Grok Compute machines as needed.

## Monitoring

Keep your hand on the platform's pulse. [Health endpoints](../../../../develop/server-management.md#server-health) show the state of every
service. [Usage Analysis](../../../../govern/audit/usage-analysis.md) shows who does what, which functions and queries
run, which errors users hit, and how the database and server perform. **Log sync** pushes logs, alerts, and a heartbeat
to Amazon CloudWatch, Google Cloud Logging, or any OpenTelemetry collector, so your monitoring system can page you
when a service degrades.

[Learn more](../../../../govern/manage-enterprise.md#monitor-the-platform).

## Audit

Know exactly what was done, when, and by whom by auditing user activity. Optimize your processes by analyzing system
usage data.

## Integration

Datagrok is extensible by design. You can develop server-side or client-side
extensions in many languages. You can add new queryable data sources using Grok Data API,
build new viewers on top of Datagrok using Grok JS API, or develop custom applications on top of it.

## Administration

Built-in admin tools let you change hundreds of parameters and defaults that are exposed by the platform. Use jobs or
alerts to automate anything.

[Learn how to manage an enterprise instance](../../../../govern/manage-enterprise.md).
