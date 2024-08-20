---
title: "Versioning policy"
position: 4 # float position is supported
---

We follow one of the most common versioning strategies for our
products ([packages](https://datagrok.ai/help/develop/#packages), [libraries](https://github.com/datagrok-ai/public/tree/master/libraries),
[JS API](https://datagrok.ai/help/develop/js-api), [Docker images](https://hub.docker.com/u/datagrok)):
[Semantic versioning](https://semver.org/). This convention helps to

- Understand the severity of changes in each new distribution
- Avoid "[dependency hell](https://en.wikipedia.org/wiki/Dependency_hell)."

Semantic versioning has a defined structure and rules that dictate how the version changes. It consists of
`MAJOR.MINOR.PATCH`

1. **PATCH** version increments when you make backward-compatible bug fixes
2. **MINOR** version increments when you add functionality in a backward-compatible manner
3. **MAJOR** version increments when you make incompatible changes

The Semantic versioning scheme is fully documented on the [official website](https://semver.org/); see [our little cheat
sheet](#versioning-rules) to help you understand our products' versioning basics.

The software in the below instruction can be
either [Datagrok JS API](https://datagrok.ai/help/develop/js-api),
[Datagrok package](https://datagrok.ai/help/develop/#packages),
[Datagrok library](https://github.com/datagrok-ai/public/tree/master/libraries),
or [Datagrok Docker image](https://hub.docker.com/u/datagrok).

## Versioning rules

1. Use the major version zero (`0.y.z`) for initial development. The versions won't be available on the Datagrok
   platform. At this stage, anything MAY change at any time. The software is not considered stable.
2. Upgrade the software to the major version first (`1.0.0`) as soon as the software is used in production. The major
   version greater than zero indicates that users can rely on the software. Follow the versioning policy described in
   the document after the first major version is released.
3. If the software (hereinafter dependent) requires a new version of any other datagrok software (hereinafter
   dependency), for example, [JS API](https://github.com/datagrok-ai/public/tree/master/js-api)
   or [datagrok library](https://github.com/datagrok-ai/public/tree/master/libraries), change the minimum required
   version of the dependency.
   If it is not yet released, set the dependency version to the local relative path of the dependency, for
   example, `"datagrok-api": "../../js-api"`. The dependent will be released once a new dependency version is released.
   The CI system will take care of the changes.

4. Use the patch version increment for bug fixes, performance improvements, refactoring, and security updates.
5. Increment minor version if new functionality is introduced in the software.
6. Avoid breaking compatibility in the libraries, JS API, and package API. If you need to break the behavior of the
   existing
   functionality: mark old method as obsolete in a minor version, check
   that everyone switched to the new method, then delete an outdated method in a new major version.
   If you want to remove libraries, JS API, or package API functionality which nobody uses: double-check with the
   colleagues
   and in the code and then release a new minor version with a removed method.
7. If breaking changes were introduced in the library, JS API, or package API, use major version increment
8. Use your common sense if you cannot find any examples neither in the [versioning rules](#versioning-rules) nor
   the [official SemVer documentation](https://semver.org/). And remember to update the
   dependencies if needed.

## Release process

Before publishing the release, the new version goes through all [quality assurance](../qa/quality-assurance.md) operations,
including the QA Engineer's manual testing of critical changes. The full process of release can be found
in [wiki](ci-flow.mdx#release-cicd-flow).

How to request the release:

1. **PATCH** release can be requested by anyone on the team. The main difference between other releases is that patch
   release includes only specific changes. To ask for the patch release, send the DevOps engineer the commits that need
   to be included.
2. **MINOR** release is the most common release type with the most recent bug fixes and features. We release the minor
   version every two weeks.
3. **MAJOR** release contains breaking changes. We avoid introducing incompatible changes, so it is infrequent. Only
   managers can request major releases.

## Bleeding-edge

Bleeding-edge releases (`bleeding-edge` Docker image tag) are unstable. We do not recommend them for production usage.
It may contain breaking changes, unstable features, and serious errors, which could cause crashes.

The image is built overnight and released to [the development environment](https://dev.datagrok.ai). It includes all the
latest developments, which
help develop [packages](../../develop/develop.md#packages).

## Check current version

Every Docker container has information about the version and branch it was built on:

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

## More information

- [Changelog policy](changelog-policy.md)
- [Release History](../../deploy/releases/release-history.md)
- [Quality assurance](../qa/quality-assurance.md)
- [CI/CD Flow](ci-flow.mdx)
