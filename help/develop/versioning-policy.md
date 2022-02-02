<!-- TITLE: Versioning policy -->
<!-- SUBTITLE: -->

# Versioning policy

## Docker Image Versioning policy

The Datagrok versioning policy requires to include a major, minor, and test version number and look
like: `MAJOR.MINOR.TEST`.

### Major release

A major release is made every time the code is released on [public.datagrok.ai](https://public.datagrok.ai). Before the
release, the new version goes through all [quality assurance](admin/quality-assurance.md) operations, including manual
testing of critical changes by the QA Engineer.

It may contain changes that break backward compatibility with previous versions or represent fundamental changes to
concepts.

For a major release, the `MAJOR` component must be incremented by one, and the `MINOR` and `TEST`
components must be set to zero.

### Minor release

Minor releases are the most common type of release which has the most recent features and bugfixes. Every release is
tested on the dev environment and passes all [quality assurance](admin/quality-assurance.md) requirements.

For a minor release, the `MINOR` component must be incremented by one, and the `TEST` component must be set to zero.
The `MAJOR` component must remain unchanged.

### Test release

Test releases are unstable. We do not recommend them for production usage.

It may contain breaking changes, unstable features, serious errors, could cause crashes. Also, it may not include all
the features that are already in other versions.

For a test release, the `TEST` component must be incremented by one, and the `MINOR` and `MAJOR`
components must remain unchanged.

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

See also:

* [Release History](release-history.md)
* [Stable Releases](release-stable.md)
