---
title: "Versioning policy"
---

## Docker Image Versioning policy

The Datagrok versioning policy requires to include a major, minor, and test version number and look
like: `MAJOR.MINOR.PATCH-<BUILD>`.

### Major release

A major release is made every time the code is released on [public.datagrok.ai](https://public.datagrok.ai). Before the
release, the new version goes through all [quality assurance](../quality-assurance.md) operations, including manual
testing of critical changes by the QA Engineer.

It may contain changes that break backward compatibility with previous versions or represent fundamental changes to
concepts.

For a major release, the `MAJOR` component must be incremented by one, and the `MINOR` and `TEST`
components must be set to zero.

The [CI/CD Flow for patch releases](../../advanced/ci-flow.md#major-release) includes all checks, and deployment to all
environments.

### Minor release

Minor releases are the stable release which has the most recent features. Every release is tested on the test
environment and passes all [quality assurance](../quality-assurance.md) requirements.

For a minor release, the `MINOR` component must be incremented by one, and the `PATCH` component must be set to zero.
The `MAJOR` component must remain unchanged.

The [CI/CD Flow for patch releases](../../advanced/ci-flow.md#minor-release) includes all checks, and deployment to dev
and test
environments.

### Patch release

Patch releases are the most common type of release which has the most recent bugfixes and small new features. Every
release is tested on the dev environment and passes all [quality assurance](../quality-assurance.md) requirements.

For a patch release, the `PATCH` component must be incremented by one, and the `MINOR` and `MAJOR`
components must remain unchanged.

The [CI/CD Flow for patch releases](../../advanced/ci-flow.md#patch-release) includes all checks, and deployment to dev
environment.

### Build release

Build releases are unstable. We do not recommend them for production usage.

It may contain breaking changes, unstable features, serious errors, could cause crashes. Also, it may not include all
the features that are already in other versions.

For a test release, the `BUILD` component must be added to version, and the `PATCH`, `MINOR` and `MAJOR`
components must remain unchanged.

The [CI/CD Flow for build releases](../../advanced/ci-flow.md#build-release) is simple and does not include most of the
checks.

### Version information in Docker image

Every Docker container has information about version and branch it was build on:

1. In Docker Labels

    ```bash
    docker inspect <image>:<version> --format '{{ json .Config.Labels }}'
    {
      "branch": "<branch>",
      "commit-main": "<commit>",
      "commit-public": "<commit>"
    }
    ```

2. In Docker container file system

    ```bash
    cat /home/grok/git-info.txt
    Branch main: <branch>
    Commit main: <commit>
    Commit public: <commit>
    ```

### More information

* [Release History](release-history.md)
* [Stable Releases](release-stable.md)
* [CI/CD Flow](../../advanced/ci-flow.md)
