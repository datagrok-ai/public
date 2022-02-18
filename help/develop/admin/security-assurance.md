<!-- TITLE: Security assurance -->
<!-- SUBTITLE: -->

# Security assurance

## Architecture

Datagrok was designed with the security as a main objective, we should be in complete control of the data.

### Credentials

[Security credentials](../../govern/security.md#credentials) are used to gain access to external resources. For example,
database connections typically require a pair of login and password.

Data connection credentials are managed
using [Datagrok Credentials Management Service](../../govern/security.md#credentials). All credentials are encrypted and
stored securely in [storage](../../govern/security.md#credentials-storage)

### Authentication

Datagrok supports different types of authentication including but not limited to login-password, OAUTH, LDAP.

Enterprise customers might prefer to use a custom SSO (single sign-on) scheme. We can accommodate these needs by
developing a customer-specific integration.

More information about Datagrok authentication capabilities can be found
on [Authentication](../../govern/authentication.md) page.

### User permissions

The platform has a flexible access control mechanism that lets administrator define people (or groups of people) that
are allowed to execute actions against different entities, based on the entity attributes. For instance, it is possible
to define a group of people who would be able to open dashboards, but would not have access to the underlying
connection. See
[Authorization System](../../govern/authorization.md) for details.

### Vulnerability remediation

Datagrok platform is written in Dart, which is a language developed and supported by Google. We
follow [Dart Security Best Practises](https://dart.dev/security#best-practices) to minimize the risk of introducing a
vulnerability.

Every release is tested according to our [quality-assurance procedures](quality-assurance.md) which includes
vulnerability check by [Snyk](https://snyk.io/) and [Grype](https://github.com/anchore/grype/). All critical
vulnerabilities are remediated in 2 weeks, high vulnerabilities in 1 moth and lower-level vulnerabilities in 3 month.

Every release contains the most recent dependency packages which saves us time with vulnerability remediation and
provides the best quality experience to the user.

Datagrok packages are also tested
using [CodeQL](https://codeql.github.com/) ([results are available publicly](https://github.com/datagrok-ai/public/security/code-scanning?query=tool%3ACodeQL))
and [Grype](https://github.com/anchore/grype/) ([results are available publicly](https://github.com/datagrok-ai/public/security/code-scanning?query=tool%3A%22Anchore+Container+Vulnerability+Report+%28T0%29%22))

### Client-server interactions

Datagrok client-side uses HTTP rest API to interact with server-side. Authentication token must be passed to access all
features.

## Infrastructure

Datagrok consist of Docker containers, [database](infrastructure.md#database)
and [persistent file storage](infrastructure.md#storage).

### Encryption in-transit

All client-server communications use the [HTTPS](https://en.wikipedia.org/wiki/HTTPS) protocol, which means it is secure
and encrypted.

### On-premise deployment

Enterprises typically prefer on-premise deployment for multiple reasons, such as security, ability to easily access
internal data, and other features such as integration with the enterprise
[authentication](../../govern/authentication.md) methods. Regarding Datagrok infrastructure it can be easily done. In
case of in-house deployment we rely on internal company policies.

#### CloudFormation deployment

To simplify deployment with all security policies taken into consideration we created
an [ECS CloudFormation deployment template](deploy-amazon-cloudformation.md).

CloudFormation Template is tested by [Snyk](https://snyk.io/) on every change. The results
are [available publicly](https://github.com/datagrok-ai/public/security/code-scanning?query=tool%3A%22Snyk+IaC%22)

In the resulted stand all communications are restricted
by [Security Groups](https://docs.aws.amazon.com/vpc/latest/userguide/VPC_SecurityGroups.html), services can not be
accessed directly, all requests goes through Application Load Balancer.

#### Encryption at-rest

For AWS deployment, we rely on Amazon's built-in encryption for
[RDS](https://docs.aws.amazon.com/AmazonRDS/latest/UserGuide/Overview.Encryption.html)
, [S3 buckets](https://docs.aws.amazon.com/AmazonS3/latest/dev/bucket-encryption.html),
and [server-side encryption of ECS ephemeral storage](https://aws.amazon.com/blogs/containers/introducing-server-side-encryption-ephemeral-storage-using-aws-fargate-managed-keys-aws-fargate-platform-version-1-4/)
.

## More information:

* [Infrastructure](infrastructure.md)
* [CloudFormation deployment](deploy-amazon-cloudformation.md)
* [Security](../../govern/security.md)
