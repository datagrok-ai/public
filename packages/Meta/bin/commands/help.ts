const HELP = `
Usage: grok-meta <COMMAND>

Datagrok's Meta package management tool

Commands:
    add       Add a package change to Meta package

To get help on a particular command, use:
    grok <COMMAND> --help
`

const HELP_ADD = `
Usage: grok add --package <name> --ver <semver> --repository <repositoryJson> [--description <string>] [--dep <dependency_name> --depver <semver>] [--category <category>]

--package     Required  Set package name to add/change in Meta package
--ver         Required  Set package version to add/change in Meta package
--repository  Required  Set package repository in Meta package
--description Optional  Set package description to add/change in Meta package
--dep         Optional  Set package dependency to add/change in Meta package
--depver      Optional  Set package dependency version to add to Meta package
--category    Optional  Set package category in Meta package
`;

export const help = {
  add: HELP_ADD,
  help: HELP
};
