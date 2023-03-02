---
title: "Quality assurance"
---

Datagrok is an incredibly powerful platform. To give users the best possible
performance, we have created a number of unique, proprietary technologies, such
as [in-memory columnar database](architecture.md#data-engine) and fast
[visualizations](architecture.md#viewers). To provide the seamless
end-to-end-experience, we took ownership of [data access](../../home.md#access),
[data governance](../../home.md#access),
[data exploration](../../home.md#explore),
[scientific computing](../../compute/scripting.md),[machine learning and artificial intelligence](../../home.md#explore),
and [collaboration](../../home.md#share).
To make sure our enterprise customers can work with their data in a secure and
efficient manner, we have built features like
[authentication](../../govern/authentication.md),
[audit](../../govern/audit.md), and
[role management system](../../govern/authentication.md). The list goes on and on.

This sort of power comes at the cost of significantly increasing the platform's
complexity. And our QA processes, along with the CI/СD concept, maintain the
reliability and functionality of platform core and external systems that we
integrate with dependencies and interactions they bring.

## Continuous integration and deployment system

To ensure quick and reliable code delivery, we apply the
[CI/CD flow](../advanced/ci-flow.md)
using the following tools:

* We use Git as a VCS, and we follow the
[Versioning policy](../admin/releases/versioning-policy.md) for Datagrok
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
  5. [Snyk](https://snyk.io/) and [Grype](https://github.com/anchore/grype/)
  check for vulnerabilities
* Docker Hub stores built artifacts, such as compiled code and packages, that
are generated during the build process.
* Ansible automates the process of provisioning and configuring servers for
deployment.
* We use Docker to run platform versions in isolated environments for
[automated](#automated-testing) and [manual](#manual-testing) testing while
deploying releases. Every time a docker image is built, the following tasks get
executed:
  1. Deploying the server on dev environment
  1. Integration tests
  1. API and UI Tests
  1. Package Tests
  1. Performance benchmarks
  1. Stress testing
* Docker Swarm, AWS ECS, and K8S are our orchestration tools for the deployment
and scaling of containers.
* We use Prometheus to track the performance and behavior of applications in
production.

## Automated testing

We use automated testing to identify issues more quickly and ensure the software
meets expected quality standards. With the help of CI/CD, automated testing is
performed on code changes as they are integrated, making it possible to identify
bugs and other issues in the software's functionality. All the results of the
automated tests are kept in the test tracking system. There are several types of
automated tests that we use:

* **Unit tests**. Our core libraries contain over 2,000 unit tests that execute
automatically in the continuous integration process, enabling us to detect
regression bugs.* **Integration tests**. These tests cover interactions
between different modules, while unit tests tend to address isolated modules.
* **JS API tests**. We have a  a regular
[Datagrok package](https://github.com/datagrok-ai/public/tree/master/packages/ApiTests)
containing a suite of tests for the [JS API](../js-api.md). These tests are
executed automatically on a server using a headless browser mode.
* **UI tests**. We have over 100 complex UI tests for automated testing of the
platform.  Each test is associated with a separate story, which describes the
intended objective, actions, and expected outcome of the user. To emulate user
input, we use Selenium. It also checks the expected result.
* [**Package tests**](../how-to/add-package-tests.md). These tests include
[unit tests](../how-to/add-package-tests.md/#adding-unit-tests) and
[function tests](../how-to/add-package-tests.md/#testing-functions). We provide
function tests for packages, scripts, and APIs, as they all utilize the
concept of
[functions](../../datagrok/functions/functions.md).
* **Performance benchmarks**. We run performance benchmarks as part of the build
and keep track of the results to ensure that the platform only gets better
with time.
* **Stress tests**. To ensure that the platform is stable under heavy load, we
perform automated stress testing regularly. The real workflow is emulated by
spawning many virtual machines that have Selenium installed on them and
executing selected UI tests on them. This helps us identify steps responsible
for the performance degradation.  For more details on different environments,
tests, and interpretation of results, see
[stress testing results](stress-testing-results.md).

## Manual testing

We use smoke tests and structured tests for complex scenarios in two primary
cases:

* **Testing new or fixed functionality**. Whenever new functionality is
developed or existing functionality is modified, the QA engineer executes test
cases that validate the expected behavior of the software. This testing
ensures that the functionality is working as expected and any defects found
are logged and fixed before the release of the software.
* **Release testing** Before the software is released, the QA engineer executes
a comprehensive set of test cases to ensure the software is ready for release
and meets the specified quality standards.

## QA monitoring tools

### Test tracking system

In order to monitor the results of test execution, we have developed a unique
application using the Datagrok platform. This application provides a visual
representation of the application's current state and enables QA engineers
promptly access, execute, and report the completion status of stories.
Additionally, it lets us quickly navigate to the corresponding JIRA issues as we
have integrated JIRA with Datagrok's OpenAPI capabilities.

![Test Tracking System](test-tracking-system.png)

### Test manager

Is a tool within the Datagrok platform that provides a convenient interface for
selecting and running packages unit tests with additional capabilities for
exploring the test results. **Test manager** is a component of the
[**DevTools** package](https://github.com/datagrok-ai/public/tree/master/packages/DevTools).
See [Test Manager](../how-to/test-packages.md#test-manager) to learn more.

### Usage Analysis package

It is the
[Datagrok package](https://github.com/datagrok-ai/public/tree/master/packages/UsageAnalysis)
designed to provide insight into how the Datagrok platform is utilized. Usage
Analysis offers various viewers, enabling diverse information pertaining to
usage statistics, platform events, and errors.
See [Usage Analysis](../../govern/usage-analysis.md)
to learn more.

### Grok Tester

## Issue tracking

Issue tracking is an essential component of any software development process,
and we work on a reliable process to identify and address issues promptly. At
Datagrok, we have multiple sources of getting feedback, and we track issues we
achieve from our clients through emails, Teams, meetings, and clients' QA.
Additionally, we use manual tests, automatic tests, monitoring, a logging
system, and our [community](https://community.datagrok.ai/) to gather feedback
and identify issues.

Users can report an issue directly in one of our issue
tracking systems. We use [JIRA](https://reddata.atlassian.net/) for internal
issues, integrated with BitBucket, and  [GitHub
Tracker](https://github.com/datagrok-ai/public/issues) for all externally
reported issues. We recommend reporting issues in JIRA if you have access to it,
but if you don't, use GitHub Tracker. 

Tracking the issues, we stay on
top of the platform's current status, including critical issues, promised
features, and upcoming changes. You can also track our latest updates in [Release notes](../admin/releases/release-history.md).

See also:

* [Architecture](architecture.md)
* [Infrastructure](infrastructure.md)
