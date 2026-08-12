// This package ships no client code. It exists so the platform manages the five kernel
// gateways as on-demand Docker containers — see dockerfiles/<language>/container.json.
//
// The file is here because the platform requires an entry point for every package:
// _publishPackage rejects a package with neither package.js nor webpack.config.js, so
// without it the package downloads, fails to publish, and is then skipped permanently
// because the install query filters on `packages.error is null`.
//
// Deliberately not a webpack package: declaring webpack.config.js would instead require a
// built dist/package.js, which there is equally nothing to build.
