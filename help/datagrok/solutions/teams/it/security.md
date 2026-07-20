---
title: "Security"
description: Overview of Datagrok's security architecture, authentication, credential management, and vulnerability remediation.
keywords:
  - authentication and authorization
  - credential management
  - vulnerability scanning
  - cve remediation
  - openvex
  - encryption at rest
  - on-premises security
  - codeql scanning
---

Datagrok was designed with security as one of the primary objectives. The platform gives users complete control of their
data, letting them decide exactly what is shared with whom, and using the best security practices across the solution.
At the same time, it empowers customers' IT to manage the platform.

Most Datagrok customers are pharma and biotech companies that deploy on-premises
because the platform works best when integrated with internal databases and systems,
which often contain proprietary chemical compounds. On-premises deployment requires
conforming to internal IT regulations but simplifies security assessment because
the platform runs in a controlled, IT-managed environment.

Datagrok's in-browser data processing improves security. Unlike other products that
require uploading datasets to a server before you can work with them, Datagrok can
perform many operations directly in the browser, reducing complexity and attack
surface. When server-side processing is needed, the platform handles it
transparently. Plugins can also run algorithms against data already loaded in the
browser, keeping sensitive data local.

[Apps](../../../../develop/develop.md#packages) built on the platform are more secure than traditional web applications
because they reuse the platform's security system instead of implementing their own.

## Authentication

Datagrok supports different types of authentication including but not limited to login-password, OAUTH, and LDAP.

Enterprise customers might prefer to use a custom SSO (single sign-on) scheme. We can accommodate these needs by
developing a customer-specific integration.

More information about Datagrok authentication capabilities can be found on
the [Authentication](../../../../govern/access-control/access-control.md#authentication) page.

## Credentials

[Security credentials](../../../../govern/access-control/access-control.md#credentials-management-system) are used to gain access to external resources. For example,
database connections typically requires a pair of login and password.

Data connection credentials are managed
using [Datagrok Credentials Management Service](../../../../govern/access-control/access-control.md#credentials-management-system). All credentials
are [encrypted and stored securely](../../../../govern/access-control/access-control.md#credentials-storage).

In case of AWS deployment, you can bypass Datagrok Credential Management Service and use
[AWS Secret Manager](../../../../govern/access-control/data-connection-credentials.md#aws-connection) instead.

Once a connection is set up, access to it (either `use` or `edit`) is subject to [user permissions](#user-permissions).

## User permissions

The platform has a flexible access control mechanism that lets administrators define users (or groups of users) that are
allowed to execute actions against different entities, based on the entity attributes. For instance, it is possible to
define a group of people who would be able to open dashboards, but would not have access to the underlying connection.
See
[Authorization System](../../../../govern/access-control/access-control.md#authorization) for details.

## Vulnerability remediation

Datagrok server core is written in Dart, a language developed and supported by Google. We
follow [Dart Security Best Practices](https://dart.dev/security#best-practices) to minimize the
risk of introducing a vulnerability. Every release ships the most recent dependency packages,
which reduces remediation effort and gives you the best-quality experience.

We continuously scan everything we ship for known vulnerabilities (CVEs), as part of our
[quality-assurance procedures](../../../../develop/qa/quality-assurance.md):

* **Docker images** — every shipped image (core services, package plugins, and base/tools
  images) is scanned with
  [Google Cloud Artifact Analysis](https://cloud.google.com/artifact-analysis/docs/container-scanning-overview),
  both weekly and when the image is built or published. Results are published openly as VEX (see below).
* **NPM packages** — the production dependency tree of every published Datagrok npm package
  (plugins, libraries, `datagrok-api`, `datagrok-tools`) is audited with
  [npm audit](https://docs.npmjs.com/cli/commands/npm-audit), both weekly and at publish time.
  Results are published openly as VEX (see below).
* **Source code** is scanned with [CodeQL](https://codeql.github.com/).
* **Dependencies and packages** are scanned with [Grype](https://github.com/anchore/grype/)
  ([results](https://github.com/datagrok-ai/public/actions/workflows/security_scan_anchore.yaml)).
* **Deployment templates** — the
  [CloudFormation template](../../../../deploy/aws/deploy-amazon-eks.mdx) is checked by
  [Snyk](https://snyk.io/) ([results](https://github.com/datagrok-ai/public/actions/workflows/iaac.yaml)).

We remediate critical vulnerabilities within two weeks, high within one month, and
lower-severity ones within three months.

### Vulnerability disclosures (VEX)

We publish the scan results openly — Docker images and npm packages alike — so you can assess
Datagrok against your own policies. Findings are provided as machine-readable
[VEX](https://www.cisa.gov/resources-tools/resources/minimum-requirements-vulnerability-exploitability-exchange-vex)
([OpenVEX](https://openvex.dev/)) documents, with CSV and human-readable HTML reports alongside.
Browse them from the index:

* **Index** — [data.datagrok.ai/vex/index.html](https://data.datagrok.ai/vex/index.html)
  (machine-readable: [index.json](https://data.datagrok.ai/vex/index.json)). Lists every image
  and npm package with its current version, vulnerability counts by severity, and links to its
  reports.
* **Composite VEX** — one document per group:
  [core services](https://data.datagrok.ai/vex/core.json),
  [package images](https://data.datagrok.ai/vex/packages.json),
  [tools and base images](https://data.datagrok.ai/vex/tools.json), and
  [npm packages](https://data.datagrok.ai/vex/npm.json).
* **Per image** — `https://data.datagrok.ai/vex/<image>/latest.json` (OpenVEX), with
  `latest.csv` and `latest.html` alongside each.
* **Per npm package** — `https://data.datagrok.ai/vex/npm/<name>/latest.json`, where `<name>`
  is the package name without the `@`, with `/` replaced by `-`
  (for example, `@datagrok/chem` → `npm/datagrok-chem`). The audit covers production
  dependencies only: `devDependencies` never reach a consumer install.

## Client-server interactions

Datagrok client-side uses HTTP rest API to interact with server-side. Authentication token is passed with each call, and
the server checks its validity before proceeding.

## Infrastructure

Datagrok consist of Docker containers, [database](../../../../develop/under-the-hood/infrastructure.md#1-core-components)
and [persistent file storage](../../../../develop/under-the-hood/infrastructure.md#1-core-components).

### Encryption in-transit

All client-server communications use the [HTTPS](https://en.wikipedia.org/wiki/HTTPS) protocol, which means it is secure
and encrypted.

### Encryption at-rest

For AWS deployment, we rely on Amazon's built-in encryption for
[RDS](https://docs.aws.amazon.com/AmazonRDS/latest/UserGuide/Overview.Encryption.html)
, [S3 buckets](https://docs.aws.amazon.com/AmazonS3/latest/dev/bucket-encryption.html),
and [server-side encryption of ECS ephemeral storage](https://aws.amazon.com/blogs/containers/introducing-server-side-encryption-ephemeral-storage-using-aws-fargate-managed-keys-aws-fargate-platform-version-1-4/)
.

### On-premise deployment

Enterprises typically prefer on-premise deployment for multiple reasons, such as security, ability to easily access
internal data, and other features such as integration with the enterprise
[authentication](../../../../govern/access-control/access-control.md#authentication) methods. In case of on-premise deployment, we rely on the internal
company policies.

#### CloudFormation deployment

To simplify deployment with all security policies taken into consideration, we created
an [EKS CloudFormation deployment template](../../../../deploy/aws/deploy-amazon-eks.mdx).

CloudFormation Template is tested by [Snyk](https://snyk.io/) on every change. The results
are [available publicly](https://github.com/datagrok-ai/public/actions/workflows/iaac.yaml).

In the resulting deployment, [Security Groups](https://docs.aws.amazon.com/vpc/latest/userguide/VPC_SecurityGroups.html)
restrict all communications. Services cannot be accessed directly — all requests
go through the Application Load Balancer.

## More information

* [Architecture](../../../../develop/under-the-hood/architecture.md)
* [Infrastructure](../../../../develop/under-the-hood/infrastructure.md)
* [CloudFormation deployment](../../../../deploy/aws/deploy-amazon-eks.mdx)
* [Security](../../../../govern/access-control/access-control.md#credentials-management-system)
