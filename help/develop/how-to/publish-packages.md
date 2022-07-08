<!-- TITLE: Publish packages -->

# Publish package

To make your package publicly available, you need to publish a package to NPM. Then you will be able to install it to
any Datagrok instance using the package manager.

## Public packages

Commit and push a new version in `package.json` for a package to
the [repository master branch](../../collaborate/public-repository.md). Then GitHub Actions will handle the rest.

1) Clone public [repository](../../collaborate/public-repository.md)
2) Make changes in the package or create a new package in the `packages/` directory
3) [Test your changes locally](test-packages.md#local-testing)
4) Bump version for a package in `package.json`. The package version must comply with Semantic Versioning. All packages
   in the [repository](../../collaborate/public-repository.md) must be located in `@datagrok` scope.
5) Commit and push all changes, including `package.json`, to the [repository](../../collaborate/public-repository.md)
6) On the push action GitHub Actions will run workflow:
    1) Build the package first
    2) Then, execute unit tests for the package before publishing.
    3) After successful build and test, the package is pushed to [NPM registry](https://www.npmjs.com/)
7) Check that the GitHub Actions workflow finished successfully
   ![GitHub Action publish status](github-actions-publish-status.png)
8) Now, you can install a new package version to the Datagrok platform using the package manager

### Trigger GitHub Actions manually

If an error occurred for the action triggered by the commit, it is possible to trigger the action manually. GitHub will
check the package version from `package.json` in the NPM repository. If the version is absent, it will run the
publishing job. Otherwise, it will run the test and build jobs instead.

1) Use [Packages workflow](https://github.com/datagrok-ai/public/actions/workflows/packages.yml)
2) Press `run workflow` and set packages list to publish separated with spaces, for example: `Demo Tutorials`. Use
   the `master` branch. Then `Run workflow`
3) Check that the GitHub Actions workflow finished successfully

The same steps can be applied
for [Libraries workflow](https://github.com/datagrok-ai/public/actions/workflows/libraries.yaml)
, [datagrok-tools package](https://github.com/datagrok-ai/public/actions/workflows/tools.yml)
and [datagrok-api package](https://github.com/datagrok-ai/public/actions/workflows/js-api.yml)

## Private packages

If the package can not be pushed to the NPM registry, you can publish it to the platform directly.

First, you need to [configure local environment](set-up-environment.md)  to access your Datagrok instance. Then, you can
publish the package to the platform. See [this guide](../develop.md#publishing) for details.

```shell
grok publish <HOST>
```

If you want to use this package with other users, share the package with the [group](../../govern/group.md).

![Share package](share-package.png)
