<!-- TITLE: Add package tests -->

# How to add package tests

Packages developed for the platform should be tested. Datagrok supports several
mechanisms for testing purposes. This article provides instructions for
utilizing these mechanisms in the package code.

## Adding unit tests

We use a custom test framework that is similar to `Jest`, and it supports
Datagrok API. To learn more about that, see our [guide for
packages](https://github.com/datagrok-ai/public/blob/master/packages/GUIDE.MD#tests).

To add tests and `Jest` configuration for local testing, run the following
command:

```shell
grok create <package-name> --jest
```

If you need to add tests support to an existing package, run this command
instead:

```shell
cd <package-name>
grok add tests
```

Once you get package dependencies (`npm install`), you can start writing tests.
First, create a folder for your test files:

```shell
cd <package-name>/src
mkdir tests
```

Each file typically includes tests for one category. Here is an example:

```js
import { category, expect, test } from '@datagrok-libraries/utils/src/test';


category('Examples', () => {
  test('Success', async () => {
    expect(1, 1);
  });

  test('Fail', () => {
    throw 'Exception';
  });
});
```

Next, make sure to import your test files in `src/package-test.ts`. After that,
build and publish your package. There are several options to runs the tests:
[locally using Jest](test-packages.md#local-testing), via [DG
console](test-packages.md#running-tests-in-the-console), and via [Test
manager](test-packages.md#test-manager). All public packages in the
[repository](../../collaborate/public-repository.md) are tested using GitHub
Actions on every commit. There is an option to trigger GitHub Actions
[manually](test-packages.md#trigger-github-actions-manually), if something goes
wrong during the auto-check.

For some real-life examples, please refer to the
[Chem](https://github.com/datagrok-ai/public/tree/master/packages/Chem) package,
which has `Jest` properly configured and all tests written properly.

## Skipping tests

If a test fails for some reason, you can skip it using the skipReason parameter
(specify a reason for skipping the test, for example,
the associated Jira issue key or GitHub issue number):

```js
test('Skipped', async () => {
  expect(1, 11);
}, {skipReason: 'GROK-99999'});
```

## Testing functions

Every package utilizes the concept of [functions](../../datagrok/functions/function.md).
Tests cases can be added directly to a function's annotation. Afterwards, the metadata
is used to test package functions automatically. Use the `test` parameter to add
test cases. A test is essentially any [grok script](../../datagrok/grok-script.md)
expression that will be evaluated to a boolean value. Here are some examples:

```ts
//name: square
//input: int x
//output: int y
//test: square(1) == 1
//test: square(2) == 4
//test: square(3) == 9
export function square(x: number): number {
  return x ** 2;
}
```

Functions with the following input/output parameter types can be tested using
the `test` annotation:

| Parameter type | Support | Input/output example                         |
|----------------|---------|----------------------------------------------|
| int            | &check; | `test: f(123) == 246`                        |
| double         | &check; | `test: 10 < f(12.5) && f(12.5) < 20`         |
| bool           | &check; | `test: f(true)`,<br>`test: returnsBool()`    |
| string         | &check; | `test: f("a") == "b"` (use `""` for strings) |
| datetime       | &check; | `test: f("1/1/2020") == Date(2020, 1, 1)`,<br>`test: DateDiff(DateTime(2020, 1, 1, 3, 0, 0, 0), "2020-01-01") == 10800000` |
| map            | &check; | `f({"a": 10, "b": "c"})`,<br>`test: f("abc").testparam == "abc"` |
| dataframe      | * | Works via an additional function call<br>`f(getDataframe())` |
| column_list    | * | Requires a dataframe parameter<br>`f(getDataframe(), ["colname1", "colname2"])`|
| column         | * | Works via an additional function call |
| file           | * | Works via an additional function call |
| blob           | * | Works via an additional function call |

See also:

- [Packages](../develop.md#packages)
- [Package testing](test-packages.md)
- [Grok script](../../datagrok/grok-script.md)
- [Instructions for datagrok-tools](https://github.com/datagrok-ai/public/tree/master/tools#datagrok-tools)
