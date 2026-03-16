---
name: publish-package-to-npm
description: Publish a Datagrok package to NPM (public) or directly to a platform instance (private)
---

# Publish Package

Help the user publish a Datagrok package, either publicly via NPM + GitHub Actions or privately to a specific instance.

## Usage
```
/publish-package [package-name] [--private] [--host <server>]
```

## Instructions

### Public packages (NPM via GitHub Actions)

This is the standard path for packages in the public Datagrok repository.

1. **Make changes** in the package under `packages/` in the [public repository](https://github.com/datagrok-ai/public).

2. **Test locally** before publishing:
   ```shell
   cd packages/MyPackage
   npm run test
   ```

3. **Bump version** in `package.json`. The version must follow Semantic Versioning (`X.Y.Z` or `X.Y.Z-rc` or `X.Y.Z-rc.N`). All public packages must be in the `@datagrok` scope.

4. **Commit and push** all changes including `package.json` to the master branch.

5. **GitHub Actions workflow** runs automatically on push:
   - Builds the package
   - Runs unit tests
   - On success, publishes to NPM registry

6. **Verify** the workflow succeeded in the GitHub Actions tab.

7. **Install** the new version on any Datagrok instance via the package manager.

### Triggering GitHub Actions manually

If the automated action failed, re-trigger manually:

1. Go to the [Packages workflow](https://github.com/datagrok-ai/public/actions/workflows/packages.yml).
2. Click "Run workflow", set the packages list (space-separated names, e.g., `Demo Tutorials`), select the `master` branch.
3. GitHub checks the version against NPM; if absent, it publishes; otherwise, it runs build + test only.

Other workflows available:
- [Libraries workflow](https://github.com/datagrok-ai/public/actions/workflows/libraries.yaml)
- [datagrok-tools](https://github.com/datagrok-ai/public/actions/workflows/tools.yml)
- [datagrok-api](https://github.com/datagrok-ai/public/actions/workflows/js-api.yml)

### Private packages (direct to platform)

For packages that cannot be published to NPM:

1. **Configure environment** with server credentials in `~/.grok/config.yaml`:
   ```yaml
   default: 'dev'
   servers:
     dev:
       url: 'https://dev.datagrok.ai/api'
       key: 'developer-key-here'
   ```

2. **Build and publish**:
   ```shell
   cd packages/MyPackage
   webpack
   grok publish <HOST>
   ```
   Replace `<HOST>` with the configured server name (e.g., `dev`, `local`).

3. **Share** the package with a user group if others need access. This is done in the Datagrok UI by right-clicking the package and selecting Share.

### Publish modes

| Command                        | Description                              |
|--------------------------------|------------------------------------------|
| `grok publish`                 | Debug mode to default server             |
| `grok publish dev`             | Debug mode to 'dev' server               |
| `grok publish --release`       | Public release to default server         |
| `grok publish --build`         | Run webpack before publishing            |
| `grok publish --rebuild`       | Force rebuild before publishing          |
| `grok publish --skip-check`    | Skip validation                          |

## Behavior

- Ask the user whether they are publishing a public or private package.
- For public packages, verify the version in `package.json` has been bumped before committing.
- For private packages, help configure server credentials if needed.
- Warn if the package version already exists on NPM (for public packages).
- Remind the user to test locally before publishing.
- If the GitHub Actions workflow fails, suggest manual re-trigger with specific instructions.
- Follow project coding conventions when modifying any package files.
