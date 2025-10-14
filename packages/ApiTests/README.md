# API Tests

API Tests is a [package](https://datagrok.ai/help/develop/#packages) for the [Datagrok](https://datagrok.ai) platform.

To add tests for Datagrok's JS API:

1. Create a folder under `src/` for your test files (each file typically includes tests for one category)
2. Get the package dependencies (`npm install`)
3. Import the required utilities (refer to other tests as an example):

   ```js
   import { category, expect, test } from '@datagrok-libraries/utils/src/test';
   ```

4. Write some tests
5. Import your test files in `src/package-test.ts`
6. Publish the package
7. Open Datagrok and start the console (`~` or `Windows | Console`)
8. Launch tests for a category via `ApiTests:test(category="category-name")`, e.g., `ApiTests:test(category="Layouts")`, or a specific test via `ApiTests:test(category="category-name", test="test-name")`, e.g., `ApiTests:test(category="Layouts", test="ViewLayout.toJson()")`, and wait for the results
