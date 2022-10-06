const Configuration = {
  /*
   * Resolve and load @commitlint/config-conventional from node_modules.
   * Referenced packages must be installed
   */
  extends: ['@commitlint/config-conventional'],
  /*
   * Resolve and load conventional-changelog-atom from node_modules.
   * Referenced packages must be installed
   */
  parserPreset: {
    parserOpts: {
      headerPattern: /^([\w ]*): (?:([\w]*))?(?:\[((?:[\w]+-[0-9]+)|(?:#[0-9]+))\])?!?: (.*)$/,
      headerCorrespondence: ['scope', 'type', 'ticket', 'subject'],
    }
  },
  /*
   * Resolve and load @commitlint/format from node_modules.
   * Referenced package must be installed
   */
  formatter: '@commitlint/format',
  /*
   * Any rules defined here will override rules from @commitlint/config-conventional
   */
  rules: {
    "header-max-length": [1, "always", 70],
    "body-full-stop": [1, "never", '.'],
    "body-leading-blank": [2, "always"],
    "body-empty": [1, "never"],
    "body-max-length": [0, "always", "Infinity"],
    "body-max-line-length": [1, "always", 70],
    "body-min-length": [0, "always", 0],
    "body-case": [1, "always", "sentence-case"],
    "footer-leading-blank": [2, "always"],
    "footer-empty": [1, "never"],
    "footer-max-length": [0, "always", "Infinity"],
    "footer-max-line-length": [1, "always", 70],
    "footer-min-length": [0, "always", 0],
    "header-case": [1, "always", "sentence-case"],
    "header-full-stop": [1, "never", "."],
    "header-max-length": [2, "always", 80],
    "header-min-length": [2, "always", 1],
    "references-empty": [1, "never"],
    "scope-enum": [0, "always", []],
    "scope-case": [0, "always", "pascal-case"],
    "scope-empty": [2, "never"],
    "scope-max-length": [0, "always", 25],
    "scope-min-length": [2, "always", 1],
    "subject-case": [2, "always", ['sentence-case', 'start-case', 'lower-case']],
    "subject-empty": [2, "never"],
    "subject-full-stop": [2, "never", "."],
    "subject-max-length": [2, "always", 60],
    "subject-min-length": [2, "always", 1],
    "subject-exclamation-mark": [0, "never"],
    "type-enum": [2, "always", ["Fix", "Feat", "Perf", "Refactor", "Docs", "Test", "Build", "CI", "Revert"]],
    "type-case": [0, "always", "pascal-case"],
    "type-empty": [0, "never"],
    "type-max-length": [1, "always", 10],
    "type-min-length": [0, "always", 0],
    "signed-off-by": [0, "always", "Signed-off-by:"],
    "trailer-exists": [0, "always", "Signed-off-by:"],
  },
  /*
   * Functions that return true if commitlint should ignore the given message.
   */
  ignores: [
    (commit) => commit === '',
    (commit) => commit.includes('WIP')
  ],
  /*
   * Whether commitlint uses the default ignore rules.
   */
  defaultIgnores: true,
  /*
   * Custom URL to show upon failure
   */
  helpUrl:
    'https://datagrok.ai/help/develop/admin/releases/versioning-policy',
};

module.exports = Configuration;
