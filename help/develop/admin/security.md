<!-- TITLE: Security -->
<!-- SUBTITLE: -->

# Security

Datagrok was designed with security as one of the primary objectives. The platform gives users complete control of their
data, letting them decide exactly what is shared with whom, and using the best security practices across the solution.
At the same time, it empowers customers' IT to manage the platform.

At the moment, most of the Datagrok customers are big pharma and biotech companies. Due to the nature of the platform (
you get maximum value when it is integrated with multiple databases and internal systems), and the data it deals with (
internal databases, often containing proprietary chemical compounds), our clients deploy Datagrok on-premises. On one
hand, it makes the integration more complicated due to the need to integrate with the internal IT systems and adhere to
the established regulations. On the other hand, that makes security assessment easier since the platform is used in a
controlled environment using the systems chosen by the IT.

Datagrok's revolutionary in-browser data processing makes it a lot more secure compared to other products. For instance,
in other products when you work with local files you need to upload the dataset to the server first in order to do
anything meaningful, therefore introducing additional layer of complexity and the need to manage that dataset on a
server. With Datagrok, a lot could be achieved without having to use the server at all. Another mode that Datagrok
supports is when the algorithm (in a form of a plugin) goes to the dataset already located in the browser - this is also
a more secure way to work with data than the other way around.

The [apps](../develop.md) built on top of the platform end up a lot more secure than the traditional web applications,
since they reuse the existing security system and do not have to roll out their own.

## Authentication

Datagrok supports different types of authentication including but not limited to login-password, OAUTH, and LDAP.

Enterprise customers might prefer to use a custom SSO (single sign-on) scheme. We can accommodate these needs by
developing a customer-specific integration.

More information about Datagrok authentication capabilities can be found on
the [Authentication](../../govern/authentication.md) page.

## Credentials

[Security credentials](../../govern/security.md#credentials) are used to gain access to external resources. For example,
database connections typically requires a pair of login and password.

Data connection credentials are managed
using [Datagrok Credentials Management Service](../../govern/security.md#credentials). All credentials
are [encrypted and stored securely](../../govern/security.md#credentials-storage).

In case of AWS deployment, you can bypass Datagrok Credential Management Service and use
[AWS Secret Manager](../../access/data-connection-credentials.md#secrets-managers) instead.

Once a connection is set up, access to it (either `use` or `edit`) is subject to [user permissions](#user-permissions).

## User permissions

The platform has a flexible access control mechanism that lets administrators define users (or groups of users) that are
allowed to execute actions against different entities, based on the entity attributes. For instance, it is possible to
define a group of people who would be able to open dashboards, but would not have access to the underlying connection.
See
[Authorization System](../../govern/authorization.md) for details.

## Vulnerability remediation

Datagrok server core is written in Dart, which is a language developed and supported by Google. We
follow [Dart Security Best Practises](https://dart.dev/security#best-practices) to minimize the risk of introducing a
vulnerability.

Every release is tested according to our [quality-assurance procedures](quality-assurance.md) which includes
vulnerability check by [Snyk](https://snyk.io/) and [Grype](https://github.com/anchore/grype/). All critical
vulnerabilities are remediated in 2 weeks, high vulnerabilities in 1 moth and lower-level vulnerabilities in 3 month.

Every release contains the most recent dependency packages which saves us time with vulnerability remediation and
provides the best quality experience to the user.

Datagrok packages are also tested using [CodeQL](https://codeql.github.com/)
and [Grype](https://github.com/anchore/grype/) ([results are available publicly](https://github.com/datagrok-ai/public/actions/workflows/security_scan.yml))

## Client-server interactions

Datagrok client-side uses HTTP rest API to interact with server-side. Authentication token is passed with each call, and
the server checks its validity before proceeding.

## Infrastructure

Datagrok consist of Docker containers, [database](infrastructure.md#database)
and [persistent file storage](infrastructure.md#storage).

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
[authentication](../../govern/authentication.md) methods. In case of on-premise deployment, we rely on the internal
company policies.

#### CloudFormation deployment

To simplify deployment with all security policies taken into consideration, we created
an [ECS CloudFormation deployment template](deploy-amazon-cloudformation.md).

CloudFormation Template is tested by [Snyk](https://snyk.io/) on every change. The results
are [available publicly](https://github.com/datagrok-ai/public/actions/workflows/iaac.yaml)

In the resulted stand all communications are restricted
by [Security Groups](https://docs.aws.amazon.com/vpc/latest/userguide/VPC_SecurityGroups.html), services can not be
accessed directly, all requests goes through Application Load Balancer.

## More information

* [Architecture](architecture.md)
* [Infrastructure](infrastructure.md)
* [CloudFormation deployment](deploy-amazon-cloudformation.md)
* [Security](../../govern/security.md)
