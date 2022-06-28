<!-- TITLE: Contribute to Help documentation -->

# Contribute to Help documentation

All Datagrok knowledge base is collected on [Help wiki pages](https://datagrok.ai/help). The source code is located in
the help directory in the [public repository](https://github.com/datagrok-ai/public/tree/master/help). Any commit to the
help directory will trigger help pages GitHub Actions workflow.

[GitHub Actions Help workflow](https://github.com/datagrok-ai/public/actions/workflows/help.yaml) run:

* Lint checks for the full documentation using [markdownlint](https://github.com/igorshubovych/markdownlint-cli) linter.
* Links check for the full documentation using [lychee](https://github.com/lycheeverse/lychee)
  linter
* Anchors check for the full documentation using [remark](https://github.com/remarkjs/remark-validate-links)
  linter
* If the checks are finished successfully, GitHub will convert the markdown files to HTML
  using [pandoc](https://pandoc.org/) and deploy them to the server.
    * The 'Deploy to server' step contains detailed information about changes that are made on the server.
* The result HTML help pages are also available as GitHub Actions artifact: `help_html_pages`

If any error occurs during lint checks, the deployment to the server will be canceled.

Check for errors in the GitHub Actions log, and search for the 'error' keyword.

## Run linter locally

```shell
# To check full documentation
docker run -v ${PWD}:/workdir --rm -i ghcr.io/igorshubovych/markdownlint-cli:latest --config .github/linters/.markdownlint.yaml "help/**/*.md"
docker run --init -it -v $(pwd)/help:$(pwd)/help -v $(pwd)/.lycheeignore:/.lycheeignore  lycheeverse/lychee $(pwd)/help --github-token <your_personal_token> --exclude-all-private --require-https --cache -a 429
npm install -g remark-cli remark-validate-links
remark -u validate-links --frail help/

# To check and fix full documentation using mardownlint
docker run -v ${PWD}:/workdir --rm -i ghcr.io/igorshubovych/markdownlint-cli:latest --fix --config .github/linters/.markdownlint.yaml "help/**/*.md"

# To check the exact file
docker run -v ${PWD}/help:${PWD}/help -v ${PWD}/.github/linters/.markdownlint.yaml:/workdir/.markdownlint.yaml  ghcr.io/igorshubovych/markdownlint-cli:latest --config .markdownlint.yaml "<absolut_file_path>"
docker run --init -it -v $(pwd)/help:$(pwd)/help -v $(pwd)/.lycheeignore:/.lycheeignore  lycheeverse/lychee "<absolut_file_path>" --github-token <your_personal_token> --exclude-all-private --require-https --cache -a 429

# To check and fix the exact file using mardownlint
docker run -v ${PWD}/help:${PWD}/help -v ${PWD}/.github/linters/.markdownlint.yaml:/workdir/.markdownlint.yaml  ghcr.io/igorshubovych/markdownlint-cli:latest --fix --config .markdownlint.yaml "<absolut_file_path>"
```

## Trigger GitHub Actions manually

If an error occurred for the action triggered by the commit, it is possible to trigger the action manually.

1) Use [Help workflow](https://github.com/datagrok-ai/public/actions/workflows/help.yaml)
2) Press `run workflow`.Choose the target branch. Then `Run workflow`. Note that deployment to the server is executed
   for the master branch only.
3) Check that the GitHub Actions workflow is finished successfully.
