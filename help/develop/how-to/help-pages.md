---
title: "Contribute to Help documentation"
---

All Datagrok knowledge base is collected on [Help wiki pages](https://datagrok.ai/help). The source code is located in
the help directory in the [public repository](https://github.com/datagrok-ai/public/tree/master/help). Any commit to the
help directory will trigger help pages GitHub Actions workflow.

[GitHub Actions Help workflow](https://github.com/datagrok-ai/public/actions/workflows/help.yaml) run:

* Lint checks for the full documentation using [markdownlint](https://github.com/igorshubovych/markdownlint-cli) linter.
* If the checks are finished successfully, GitHub will convert the markdown files to HTML
  using [Docusaurus](https://docusaurus.io/)
* Then [Hyperlink](https://github.com/untitaker/hyperlink) checks the result HTML files.
* If the checks are finished successfully GitHub Actions deploys the documentation to the server.
  * 'Deploy to server' step contains detailed information about changes that are made on the server.
* The result HTML help pages are also available as GitHub Actions artifact: `docusaurus`

If any error occurs during lint checks, the deployment to the server will be canceled.

Check for errors in the GitHub Actions log or in GitHub Actions run summary.

## Run docusaurus locally

[Docusaurus](https://docusaurus.io/) is an optimized site generator. Docusaurus
is a project for building, deploying, and maintaining open source project
websites easily.

### Requirements

1. Install Node JS [v18.12.x](https://nodejs.org/dist/v18.12.0/)
2. Upgrade npm to `v9.x.x`:

    ```shell
    npm install -g npm@9.x.x
    ```

3. Clone the public repository:

    ```shell
    git clone https://github.com/datagrok-ai/public.git
    ```

### Local setup

1. Open the terminal and change the current directory to the directory with docusaurus code:

    ```shell
    cd <absolute_path_to_public_repo>/docusaurus
    ```

2. Install package dependencies:

    ```shell
    npm install
    ```

3. As we use docusaurus to generate the API documentation, you need to repeat these installation steps in the `js-api` folder:

    ```shell
    cd <absolute_path_to_public_repo>/js-api
    npm install
    ```

4. To run docusaurus and test the changes, run the following command:

    ```shell
    cd <absolute_path_to_public_repo>/docusaurus
    npm run start
    ```

    It will start the server and recompile the client
    on every change you make. If the command does not work use `docusaurus start`.

5. Docusaurus will be available at <http://localhost:3000>. If you don't see your changes immediately, reload the page.
6. When you think the development is over and you want to commit the code, run the following command to check that the
   project builds without errors.
    1. If the command does not work use `docusaurus build`

    ```shell
    npm run build
    ```

    Note that the build step will take some time.

7. Now you can commit and push your changes to the repository.

## Run linter locally

Change `<absolute_path_to_public_repo>` and `<absolute_path_to_help_file>` to exact locations on your local filesystem.

Change `<your_personal_github_token>` to the value of your [personal GitHub token](https://github.com/settings/tokens).
If you do not have one yet, you can
easily [create it on GitHub](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token)
.

1. Markdown lint

    ```shell
    docker pull ghcr.io/igorshubovych/markdownlint-cli:latest
    docker run -v "<absolute_path_to_public_repo>":/workdir --rm -i ghcr.io/igorshubovych/markdownlint-cli:latest --config .github/linters/.markdownlint.yaml help/**/*.md
    docker run -v "<absolute_path_to_public_repo>":/workdir --rm -i ghcr.io/igorshubovych/markdownlint-cli:latest --fix --config .github/linters/.markdownlint.yaml help/**/*.md
    docker run -v "<absolute_path_to_help_file>":"<absolute_path_to_help_file>" -v "<absolute_path_to_public_repo>/.github/linters/.markdownlint.yaml":/workdir/.markdownlint.yaml  ghcr.io/igorshubovych/markdownlint-cli:latest --config .markdownlint.yaml "<absolute_path_to_help_file>"
    docker run -v "<absolute_path_to_help_file>":"<absolute_path_to_help_file>" -v "<absolute_path_to_public_repo>/.github/linters/.markdownlint.yaml":/workdir/.markdownlint.yaml  ghcr.io/igorshubovych/markdownlint-cli:latest --fix --config .markdownlint.yaml "<absolute_path_to_help_file>"
    ```

2. Links and anchors check.

   Internal links and anchors are checked after the help package convert to HTML. To test it locally wou should convert
   markdown to HTML using [Docusaurus](https://docusaurus.io/) and then check all links
   using [hyperlink](https://github.com/untitaker/hyperlink).

   ```shell
   npm install -g @untitaker/hyperlink
   cd docusaurus
   npm install
   npm run build
   hyperlink --check-anchors build/ --sources ../help/
   ```

## Trigger GitHub Actions manually

If an error occurred for the action triggered by the commit, it is possible to trigger the action manually.

1. Use [Help workflow](https://github.com/datagrok-ai/public/actions/workflows/help.yaml)
2. Press `run workflow`. Choose the target branch. Then `Run workflow`. Note that deployment to the server is executed
   for the master branch only.
3. Check that the GitHub Actions workflow is finished successfully.

See also:

* [Writing Docs: Content types](../../collaborate/writing-docs/content-types.md)
