---
title: "Security"
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

Datagrok server core is written in Dart, which is a language developed and supported by Google. We
follow [Dart Security Best Practices](https://dart.dev/security#best-practices) to minimize the risk of introducing a
vulnerability.

Every release is tested according to our [quality-assurance procedures](../../../../develop/qa/quality-assurance.md) which includes
vulnerability check by [Snyk](https://snyk.io/) and [Grype](https://github.com/anchore/grype/). All critical
vulnerabilities are remediated in 2 weeks, high vulnerabilities in 1 moth and lower-level vulnerabilities in 3 month.

Every release contains the most recent dependency packages which saves us time with vulnerability remediation and
provides the best quality experience to the user.

Datagrok packages are also tested using [CodeQL](https://codeql.github.com/)
and [Grype](https://github.com/anchore/grype/) ([results are available publicly](https://github.com/datagrok-ai/public/actions/workflows/security_scan_anchore.yaml)).

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
an [EKS CloudFormation deployment template](../../../../deploy/aws/deploy-amazon-eks.md).

CloudFormation Template is tested by [Snyk](https://snyk.io/) on every change. The results
are [available publicly](https://github.com/datagrok-ai/public/actions/workflows/iaac.yaml).

In the resulting deployment, [Security Groups](https://docs.aws.amazon.com/vpc/latest/userguide/VPC_SecurityGroups.html)
restrict all communications. Services cannot be accessed directly — all requests
go through the Application Load Balancer.

## More information

* [Architecture](../../../../develop/under-the-hood/architecture.md)
* [Infrastructure](../../../../develop/under-the-hood/infrastructure.md)
* [CloudFormation deployment](../../../../deploy/aws/deploy-amazon-eks.md)
* [Security](../../../../govern/access-control/access-control.md#credentials-management-system)
