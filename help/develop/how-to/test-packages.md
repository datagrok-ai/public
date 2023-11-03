---
title: "Test packages"
---

Testing is an essential part of development process. You should ensure that your product works properly at each stage of
it's lifecycle. For instance, when developing a new version of a package you should perform regression testing and
confirm that new changes haven't affected previous functionality. Each package should include a bunch of unit tests
responsible for either UI or logic underneath. And Datagrok provides various capabilities and tools to conveniently run
those tests any time during development.

The process of writing and executing tests can be significantly easier with the help of Datagrok library Utils,
which contains a powerful [set of testing utilities](https://github.com/datagrok-ai/public/blob/master/libraries/utils/src/test.ts).
From basic assertions to asynchronous handling and viewer testing, these utilities cover a wide range of testing scenarios.

<b>Writing Tests</b>

Use the `test` function from Utils library to define individual test cases, assign names, and specify test logic. Tests can be categorized and executed collectively for comprehensive coverage.

<details>
<summary>Example</summary>

```typescript
category('myCategory', () => {
  test('myTest', async () => {
    //do some checks
  });
});
```

</details>

<b>Test Assertions</b>

Within the Utils library, you have access to various test assertion functions to validate and compare different types of data. Functions like `expect`, `expectFloat`, `expectTable`, `expectObject`, and `expectArray` enable you to assert equality, floating-point number comparisons, DataFrame comparisons, object comparisons (including nested objects and arrays), and array comparisons, respectively.

<details>
<summary>Example</summary>

```typescript
expect(1 != null, true);
expectFloat(0.12345, 0.123, 0.001);
expectTable(grok.data.demo.demog(10), grok.data.testData('demog', 10));
expectObject({a: 1, b: 2}, {b: 2, a: 1});
expectArray([...Array(4).keys()], [0, 1, 2, 3]);
```

</details>

<b>Handling Asynchronous Tests</b>

Simplify testing asynchronous code with the `awaitCheck` utility. It allows you to wait for a specified condition or timeout, making asynchronous test cases easier to handle. The `delay` function is also available for introducing intentional delays during testing.

<details>
<summary>Example</summary>

```typescript
const smiles = grok.data.demo.molecules(20);
const v = grok.shell.addTableView(smiles);
await awaitCheck(() => document.querySelector('canvas') !== null, 'cannot load table', 3000);
grok.shell.info('Table loaded!');
await delay(1000);
v.close();
```

</details>

<b>Testing Viewers</b>

The [testViewer](https://github.com/datagrok-ai/public/blob/master/libraries/utils/src/test.ts#L451) utility is a valuable tool for testing Datagrok viewers. It facilitates thorough testing of viewer functionality, including dataframe transformation and serialization. By providing a DataFrame as input, testViewer allows you to examine how the viewer responds to various data scenarios, handling dataframe modifications such as row selection, data filtering, and column value changes. Furthermore, the utility validates viewer serialization by saving the viewer into a layout and loading it back, ensuring that viewers can accurately preserve their state.

<details>
<summary>Example</summary>

```typescript
const smiles = grok.data.demo.molecules(100);
await testViewer('Chem Similarity Search', smiles, {detectSemanticTypes: true});
```

</details>

<b>Handling Exceptions</b>

Ensure proper exception handling by using the [expectExceptionAsync](https://github.com/datagrok-ai/public/blob/master/libraries/utils/src/test.ts#L425) function. It allows you to test that specific actions throw the expected exceptions and offers the flexibility to perform deeper inspections of the exceptions when needed.

<details>
<summary>Example</summary>

```typescript
await expectExceptionAsync(() => grok.functions.call('nonExistingFunction'));
```

</details>

You can find examples of tests by referring to the
[APITests package](https://github.com/datagrok-ai/public/tree/master/packages/ApiTests).
This package contains a collection of test cases that can serve as illustrative models for your testing needs.

## Local testing

To test packages locally before release, you can use the [datagrok toolkit](../tools/libraries.md#datagrok-toolkit).
The easiest way is to execute the following commands:

```shell
cd <package-dir>
grok test [--host]
```

It will build your package and publish it in the debug mode to a specified host and run the tests. The `host` option
should correspond to one of the server aliases provided in the configuration file (`.grok/config.yaml`). If not given,
the default host is used for testing. The building and publishing steps can be skipped with flags `--skip-build` and
`--skip-publish` correspondingly. The results are printed to the terminal. If you want to save a test run result,
add the `--csv` flag (the report will be saved in a CSV file in the package folder). To see tests execution,
pass the `--gui` flag that disables the headless browser mode. This option can help you debug your tests.
For more detailed information and additional usage instructions, please refer to the [datagrok-tools page](https://github.com/datagrok-ai/public/tree/master/tools).

If you do not have any datagrok instance run locally, you can use [docker-compose](../../deploy/docker-compose/docker-compose.mdx) to run the stand.

### Trigger GitHub Actions manually

If an error occurred for the action triggered by the commit, it is possible to trigger the action manually.

1. Use [Packages workflow](https://github.com/datagrok-ai/public/actions/workflows/packages.yml)
2. Press `run workflow` and set packages list to test separated with spaces, for example: `Demo Tutorials`. Choose the
   target branch. Then `Run workflow`. Note that publish to the NPM registry is executed for the master branch only.
3. Check that the GitHub Actions workflow finished successfully
4. The results are available in the actions output

### Troubleshooting failed tests in GitHub Actions

When a GitHub Action completes, you can access the "Artifacts" section on its main page ("Summary"), which contains a zip archive with the results of the tests for that specific action. This archive typically includes three important files that can help you diagnose and resolve the problems:

1. `test-console-output.log`: This file contains the console log records generated during   the test execution. The console output often includes valuable information such as error messages, stack traces, and debugging statements. It's a good starting point for understanding what went wrong during the test run.

2. `test-record.mp4`: This file is a video recording of the test execution. The video recording can be immensely helpful in scenarios where tests are failing intermittently or in cases where the issue is related to the user interface or visual components of your application. Watching the test execution can provide insights into the behavior of the application and identify potential areas of concern.

3. `test-report.csv`: This file is a tabular representation of the test results. It contains information about which tests passed, which failed, and also additional details such as test duration or errors. The test report can give you a quick overview of the overall test status and pinpoint specific failing tests.

## Test manager

'Test manager' is a tool within the Datagrok platform that provides a convenient interface to select and run package
unit tests with further results exploration.
'Test manager' itself is a part of the DevTools package.

To start 'Test manager' go to top menu Tools -> Dev -> Test manager

![Test manager start](test-mngr-start.png)

Application starts showing a list of all packages containing unit tests. Inside each package, tests are divided by
category. Categories support multiple nesting (subcategories should be divided by `:`). To select a test or a category,
click on it, or use keyboard.

![Tests list](test-mngr-tests-list.png)

### Running tests

There are multiple ways you can run tests:

- by right clicking on package, category, or test and selecting `Run` from context menu
- by selecting package, category, or test and pushing `Enter`
- by selecting package, category, or test and pushing `Run` on a ribbon panel
- by selecting package, category, or test and pushing `Run` on a context panel
- individual tests can be run by double click
- you can run all tests at once using `Run all` button on the ribbon
- package, category, or test can be run by putting the corresponding url into address bar of the browser. The format is
  the following `your_server_name/apps/DevTools/TestManager/package_name/category_name/test_name`
- Progress icon is shown opposite to active test/category/package, it will end up in result icon after completion. In
  case at least one test fails within category or package the fail icon will be shown.

![Running tests](running_tests.gif)

Progress bar on the bottom of the page shows the percentage of completed tests.

![Progress bar](test_manager_progress_bar.png)

### Reviewing results

Information about test results is available via tooltip or in the context panel. Selected test, category, or package to
explore results. In case category/package contain multiple tests results are shown as a grid which can be added to
workspace for further exploration.

![Test results](test_results.gif)

## Running tests in the platform console

It is possible to run tests in the platform's console. Press the tilde key `~` to open the console or enable it from the
toolbox (`Windows | Console`). To
launch tests for a category, type `PackageName:test(category="category-name")`,
e.g., `ApiTests:test(category="Layouts")`, or add a specific test as a
parameter: `PackageName:test(category="category-name", test="test-name")`, e.g.,
`ApiTests:test(category="Layouts", test="ViewLayout.toJson()")`.

## Testing public package changes

It is always a good practice to test the changes before publishing the package.

All public packages in the [repository](../../collaborate/public-repository.md)
are tested using GitHub Actions on every commit. For every changed package GitHub creates a new separate instance of
Datagrok from the latest Datagrok docker image. Then, it publishes a new version of the package to this instance. And
then, the tests are executed on it.

The results are available in the actions output.

## More information

- [How to add package tests](add-package-tests.md)
