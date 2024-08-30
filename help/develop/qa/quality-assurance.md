---
title: "Quality assurance"
---

Datagrok is an incredibly powerful platform. To give users the best possible
performance, we have created a number of unique, proprietary technologies, such
as [in-memory columnar database](../../develop/under-the-hood/architecture.md#data-engine) and fast
[visualizations](../../develop/under-the-hood/architecture.md#viewers). To provide the seamless
end-to-end-experience, we took ownership of [data access](../../datagrok/datagrok.md#access),
[data governance](../../datagrok/datagrok.md#access),
[data exploration](../../datagrok/datagrok.md#explore),
[scientific computing](../../compute/scripting/scripting.mdx),[machine learning and artificial intelligence](../../datagrok/datagrok.md#explore),
and [collaboration](../../datagrok/navigation/basic-tasks/basic-tasks.md#share).
To make sure our enterprise customers can work with their data in a secure and
efficient manner, we have built features like
[authentication](../../govern/access-control/access-control.md#authentication),
[audit](../../govern/audit/audit.md). The list goes on and on.

This sort of power comes at the cost of significantly increasing the platform's
complexity. And our QA processes, along with the CI/СD concept, maintain the
reliability and functionality of platform core and external systems that we
integrate with dependencies and interactions they bring.

## Continuous integration and deployment system

To ensure quick and reliable code delivery, we apply the
[CI/CD flow](../../develop/dev-process/ci-flow.mdx)
using the following tools:

* We use Git as a VCS, and we follow the
[versioning policy](../../develop/dev-process/versioning-policy.md) for Datagrok
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
* Docker Hub stores built artifacts, such as compiled code and packages that
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
* Docker Swarm, AWS ECS, and K8S are our orchestration tools for the deployment
and scaling of containers.
* We use Prometheus to track the performance and behavior of applications in
production.

## Automated testing

With the help of CI/CD, automated testing is
performed on code changes as they are integrated. We also run automated tests every hour on our instances and encourage developers to use them when performing local package testing. All the results of the
automated tests are kept in the test tracking system.

There are several types of automated tests that we use:

* **Unit tests**. Our core libraries contain over 2,000 unit tests that execute
automatically in the continuous integration process, enabling us to detect
regression bugs.
* **Integration tests**. These tests cover interactions
between different modules, while unit tests tend to address isolated modules.
* **JS API tests**. We have a regular
[Datagrok package](https://github.com/datagrok-ai/public/tree/master/packages/ApiTests)
containing a suite of tests for the
[JS API](../../develop/packages/js-api.md). These tests are
executed automatically on a server using a headless browser mode.
* **UI tests**. We have over 100 complex UI tests for automated testing of the
platform.  Each test is associated with a separate story, which describes the
intended objective, actions, and expected outcome of the user. To emulate user
input, we use Selenium. It also checks the expected result.
* [**Package tests**](../../develop/how-to/add-package-tests.md). These tests include
[unit tests](../../develop/how-to/add-package-tests.md#adding-unit-tests) and
[function tests](../../develop/how-to/add-package-tests.md#testing-functions). We provide
function tests for packages, scripts, and APIs, as they all utilize the
concept of
[functions](../../datagrok/concepts/functions/functions.md).
* **Performance benchmarks**. We run performance benchmarks as part of the build
and keep track of the results to ensure that the platform only gets better
with time.
* **Stress tests**. To ensure that the platform is stable under heavy load, we
perform automated stress testing regularly. The real workflow is emulated by
spawning many virtual machines that have Selenium installed on them and
executing selected UI tests on them. This helps us identify steps responsible
for the performance degradation.  For more details on different environments,
tests, and interpretation of results, see
[stress testing results](../../datagrok/solutions/enterprise/stress-testing-results.md).

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

## Issue tracking

We work on a reliable process to identify and address issues promptly. At
Datagrok, we have multiple sources of getting feedback, and we track issues that we
achieve from our clients through emails, Teams, meetings, and clients' QA.
Additionally, we use manual tests, automatic tests, monitoring and logging
systems, and our [community](https://community.datagrok.ai/) to gather feedback
and identify issues.

You can [report an issue](report-tickets.md) directly in one of our issue
tracking systems. We use [JIRA](https://reddata.atlassian.net/) for internal
issues, integrated with BitBucket, and [GitHub
Tracker](https://github.com/datagrok-ai/public/issues) for all externally
reported issues. We recommend reporting issues in JIRA if you have access to it,
but if you don't, use GitHub Tracker.

Tracking the issues, we stay on
top of the platform's current status, including critical issues, promised
features, and upcoming changes. You can also track our latest updates in [Release notes](../../deploy/releases/release-history.md).

## QA tools

### GrokTester

Is an integrated monitoring system that runs automated package tests every hour. It's also possible to run [package testing locally](../../develop/how-to/test-packages#local-testing) before publishing.

### Test tracking system

To monitor the results of test execution, we have developed a unique
application on top of the Datagrok platform. This application visualizes the application's current state and enables QA engineers
promptly access, execute, and report the completion status of stories.
Additionally, it lets us quickly navigate to the corresponding JIRA issues as we
have integrated JIRA with Datagrok's OpenAPI capabilities.

![Test Tracking System](../../datagrok/solutions/enterprise/test-tracking-system.png)

### Test Manager

Test Manager is a tool within the Datagrok platform that provides a convenient interface for
selecting and running package unit tests with additional capabilities for
exploring the test results. **Test manager** is a component of the
[**DevTools** package](https://github.com/datagrok-ai/public/tree/master/packages/DevTools).
See [Test Manager](../../develop/how-to/test-packages.md#test-manager) to learn more.

### Usage Analysis

Usage Analysis is the
[Datagrok package](https://github.com/datagrok-ai/public/tree/master/packages/UsageAnalysis) that facilitates understanding of Datagrok platform usage through different viewers, displaying data on usage statistics, errors, and platform events.
See [Usage Analysis](../../govern/audit/usage-analysis.md)
to learn more.

### Error Reporting Tool and Error Analysis Dashboard

Errors encountered by users 




See also:

* [Architecture](../../develop/under-the-hood/architecture.md)
* [Infrastructure](../../develop/under-the-hood/infrastructure.md)
