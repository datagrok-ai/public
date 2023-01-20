---
title: "Contribute to Help documentation"
---

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

Change `<absolut_path_to_public_repo>` and `<absolut_path_to_help_file>` to exact locations on your local filesystem.

Change `<your_personal_github_token>` to the value of your [personal GitHub token](https://github.com/settings/tokens).
If you do not have one yet, you can
easily [create it on GitHub](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token)
.

1. Markdown lint

    ```shell
    # Download dependencies
    docker pull pandoc/core:latest
    npm install -g @untitaker/hyperlink

    # To run lint check
    docker pull ghcr.io/igorshubovych/markdownlint-cli:latest
    docker run -v "<absolut_path_to_public_repo>":/workdir --rm -i ghcr.io/igorshubovych/markdownlint-cli:latest --config .github/linters/.markdownlint.yaml help/**/*.md
    docker run -v "<absolut_path_to_public_repo>":/workdir --rm -i ghcr.io/igorshubovych/markdownlint-cli:latest --fix --config .github/linters/.markdownlint.yaml help/**/*.md
    docker run -v "<absolut_path_to_help_file>":"<absolut_path_to_help_file>" -v "<absolut_path_to_public_repo>/.github/linters/.markdownlint.yaml":/workdir/.markdownlint.yaml  ghcr.io/igorshubovych/markdownlint-cli:latest --config .markdownlint.yaml "<absolut_path_to_help_file>"
    docker run -v "<absolut_path_to_help_file>":"<absolut_path_to_help_file>" -v "<absolut_path_to_public_repo>/.github/linters/.markdownlint.yaml":/workdir/.markdownlint.yaml  ghcr.io/igorshubovych/markdownlint-cli:latest --fix --config .markdownlint.yaml "<absolut_path_to_help_file>"

    # To run external links checker
    docker pull lycheeverse/lychee
    docker run --init -it -v "<absolut_path_to_public_repo>/help":/help -v "<absolut_path_to_public_repo>/.lycheeignore":/.lycheeignore  lycheeverse/lychee /help --github-token "<your_personal_github_token>" --exclude-all-private --require-https --cache -a 429
    docker run --init -it -v "<absolut_path_to_help_file>":"<absolut_path_to_help_file>" -v $(pwd)/.lycheeignore:/.lycheeignore  lycheeverse/lychee "<absolut_file_path>" --github-token "<your_personal_github_token>" --exclude-all-private --require-https --cache -a 429
    ```

2. External links check

    ```shell
    # To run external links checker
    docker pull lycheeverse/lychee
    docker run --init -it -v "<absolut_path_to_public_repo>/help":/help -v "<absolut_path_to_public_repo>/.lycheeignore":/.lycheeignore  lycheeverse/lychee /help --github-token "<your_personal_github_token>" --exclude-all-private --require-https --cache -a 429
    docker run --init -it -v "<absolut_path_to_help_file>":"<absolut_path_to_help_file>" -v $(pwd)/.lycheeignore:/.lycheeignore  lycheeverse/lychee "<absolut_file_path>" --github-token "<your_personal_github_token>" --exclude-all-private --require-https --cache -a 429
    ```

3. Internal links and anchors check.

   Internal links and anchors are checked after the help package convert to HTML. To test it locally wou should convert
   markdown to HTML using [pandoc](https://pandoc.org/) and then check all links
   using [hyperlink](https://github.com/untitaker/hyperlink).

## Trigger GitHub Actions manually

If an error occurred for the action triggered by the commit, it is possible to trigger the action manually.

1) Use [Help workflow](https://github.com/datagrok-ai/public/actions/workflows/help.yaml)
2) Press `run workflow`.Choose the target branch. Then `Run workflow`. Note that deployment to the server is executed
   for the master branch only.
3) Check that the GitHub Actions workflow is finished successfully.
