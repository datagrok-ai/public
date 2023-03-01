---
title: "Quality assurance"
---

Datagrok is an incredibly powerful platform. To give users the best possible
performance, we have created a number of unique, proprietary technologies such
as [in-memory columnar database](architecture.md#data-engine) and fast
[visualizations](architecture.md#viewers).
To provide the seamless
end-to-end-experience, we took ownership of
[data access](../../home.md#access),
[data governance](../../home.md#access),
[data exploration](../../home.md#explore),
[scientific computing](../../compute/scripting.md),
[machine learning and artificial intelligence](../../home.md#explore), and
[collaboration](../../home.md#share).
To make sure our enterprise customers can work with their data in a secure and
efficient manner, we have built features like
[authentication](../../govern/authentication.md),
[audit](../../govern/audit.md), and
[role management system](../../govern/authentication.md).
The list goes on and on.

This sort of power comes at the cost of significantly increasing the platform's
complexity. And our QA processes, along with the CI/СD concept, maintain the
reliability and functionality of platform core and external systems that we
integrate with dependencies and interactions they bring.

## Continuous integration and deployment system

To ensure quick and reliable code delivery, we use a continuous deployment
pipeline via the following tools:

* We use Git as a VCS, and we follow the
[Semantic Versioning](https://semver.org/) strategy for Datagrok
[JS APIs](https://github.com/datagrok-ai/public/tree/master/js-api),
[packages](https://github.com/datagrok-ai/public/tree/master/packages),
[libraries](https://github.com/datagrok-ai/public/tree/master/libraries), and
[Docker images](https://hub.docker.com/u/datagrok).
* Jenkins automates the process of building and testing code every time a change
is made to the codebase in Git. It automatically executes the following tasks:
  1. Static code analysis by [dartanalyzer](https://pub.dev/packages/analyzer)
  2. Unit tests
  3. Compilation
  4. Packaging the platform in Docker containers
  5. [Snyk](https://snyk.io/) and
  [Grype](https://github.com/anchore/grype/) check for vulnerabilities

* Docker Hub stores built artifacts, such as compiled code and packages, that
are generated during the build process.
* Ansible automates the process of provisioning and configuring servers for
deployment.
* We use Docker to run platform versions in isolated environments for
[automated](#automated-testing) and
[manual](#manual-testing)
testing while deploying releases. Every time a docker image is built, the
following tasks get executed:
  1. Deploying the server on dev environment
  2. Integration tests
  3. UI Tests
  4. Performance benchmarks
  5. Stress testing
* Docker Swarm, AWS ECS, and K8S are our orchestration tools for the deployment
and scaling of containers.
* We use Prometheus to track the performance and behavior of applications in
production.

<!--We’ve developed our internal Logging System :
+ Context system
+ User settings
+ Datlas logging
Debug flags
Group settings
Client-side settings
Audit integration -->

## Automated testing

CI/CD enables automated testing of code changes as they are integrated, making
it possible to identify issues more quickly and ensure that software meets
expected quality standards.

Results of all automated tests are kept in the
[test tracking system](#test-tracking-system).

* **Unit tests**. Our core libraries contain more than 2,000 unit tests that serve two primary
purposes. First of all, they are produced and get used in the process of the
test-driven development. Secondly, they automatically execute as part of the
continuous integration process, enabling us to detect regression bugs. We also
provide unit tests for our packages and libraries, see
[Adding package unit tests](../how-to/add-package-tests.md/#adding-unit-tests)
* **Integration tests**. Integration tests are very similar to unit tests. While unit tests tend to
address isolated modules, integration tests cover interactions between different
modules.
* **JS API tests**. ! A separate suite of tests is devoted to the [JS API](../js-api.md). These
tests are executed automatically on a server using a headless browser mode. The
test suite is written in TypeScript, and is a regular
[Datagrok package](https://github.com/datagrok-ai/public/tree/master/packages/ApiTests).
* **UI tests**. ! For automated testing of the platform, we use Selenium, which works by
emulating user input, and then checking for the expected result. We have over
100 complex UI tests, and the number is proliferating. Each one is associated
with a separate story, which is a plain text file describing the objective of
the story (i.e., "user queries a database"), actions (i.e., "user opens the data
connection toolbox, ..."), and the expected outcome (i.e., "results of the query
are opened in the workspace"). Stories can be produced by external users, in
which case our QA engineer translates them to a script executable by Selenium.
Overall, these tests cover a significant part of the platform's functionality,
and there is a good chance that a bug will be caught by that system.
* **Function tests**. ! Test cases can be added directly to a function's
annotation. As every package, script, and API utilizes the concept of
[functions](../../datagrok/functions/functions.md), we provide function tests to
all of them. See
[testing functions](https://datagrok.ai/help/develop/how-to/add-package-tests#testing-functions)
for details.
* **Performance benchmarks**. We consider performance a critical feature of
Datagrok, so we keep track of it. We run performance benchmarks as part of the
build, and keep track of the results to make sure the platform only gets better
with time.
* **Stress tests**. To make sure the platform is stable under heavy load (either
many users working simultaneously or executing CPU-intensive computations), we
perform automated stress testing regularly. The real workflow is emulated by
spawning many virtual machines that have Selenium installed on them and
executing selected UI tests on them. There are two different passes:
  * An ordinary yet intensive user activity is emulated by random-ranging tests.
  * Each of the selected tests is executed on all machines simultaneously. This
helps us identify steps responsible for the performance degradation.

  For more details on different environments, tests, and interpretation of
results, see [stress testing results](stress-testing-results.md).

## Manual testing

* **Structured manual tests**. ! Certain UI tests are so complex that it's
impractical to develop and maintain corresponding scripts. In such cases, we
still write the story for it just like in UI tests. Our QA engineer executes
structured UI tests, and results are reported to the test tracking system.

* **Unstructured manual tests**. ! In addition to the testing outlined above,
each release is manually tested by our QA engineer on a separate release
environment. The amount of testing depends on the magnitude of the release.

## QA monitoring tools

### Test tracking system

In order to keep track of test execution results, we have built a unique
application on top of the Datagrok platform. It lets us not only visualize the
current state of the application but also allows the QA engineer to efficiently
work with stories by having the ability to quickly read them, execute, and
report the completion status. Moreover, it lets us quickly navigate to the
corresponding JIRA issues (integration with JIRA was done using Datagrok's
OpenAPI capabilities).

![Test Tracking System](test-tracking-system.png)

### Test manager

It is a tool within the Datagrok platform that provides a convenient interface
to select and run package unit tests with further results exploration. **Test
manager** itself is a part of the [**DevTools** package](https://github.com/datagrok-ai/public/tree/master/packages/DevTools).
See [Test Manager](../how-to/test-packages.md/#test-manager) to learn more.

### Usage Analysis package

It is a
[Datagrok package](https://github.com/datagrok-ai/public/tree/master/packages/UsageAnalysis)
that helps to understand how exactly the platform is being used. In the multiple
views, you can find different information about usage statistics, platform
events, and errors. See
[Usage Analysis](../../govern/usage-analysis.md) to
learn more.

## Issue tracking

Sources:
Clients (emails, Teams, meetings)
Clients’ QA
Manual tests
Automatic tests
Monitoring
Logging system
[Community](https://community.datagrok.ai/)

We use two issue tracking systems:

* [JIRA](https://reddata.atlassian.net/)
for internal issues integrated with BitBucket
* [GitHub Tracker](https://github.com/datagrok-ai/public/issues)
for all externally reported issues.

If you have access to Jira, we recommend reporting there; if you don't, use
GitHub.

* Know current platform status (critical issues, promised to clients, planned
  for the sprint, etc)
* Release notes

See also:

* [Architecture](architecture.md)
* [Infrastructure](infrastructure.md)
