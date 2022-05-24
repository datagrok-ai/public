const HELP = `
Usage: grok-meta <COMMAND>

Datagrok's Meta package management tool

Commands:
    add       Add a package change to Meta package

To get help on a particular command, use:
    grok <COMMAND> --help
`

const HELP_ADD = `
Usage: grok add --package <name> --description <string> --ver <semver> --dep <dependency_name> --depver <semver>

--package      Set package name to add/change in Meta package
--description  Set package description to add/change in Meta package
--ver          Set package version to add/change in Meta package
--dep          Set package dependency to add/change in Meta package
--depver       Set package dependency version to add to Meta package
`;

export const help = {
  add: HELP_ADD,
  help: HELP
};
