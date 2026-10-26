/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/add-new-column/table-and-column-functions-known-failures.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [powerpack.dialogs.add-new-column, GROK-20931, public#4404, public#4406]
--- */
import {test} from '@playwright/test';
import '../../bindings/add-new-column.js';
import '../../bindings/enrichment.js';
import '../../bindings/home.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldContainText, shouldHaveValue, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Add New Column: what the table and column functions still get wrong", () => {
  const session = feature(test, "features/add-new-column/table-and-column-functions-known-failures.feature", import.meta.url);
  test("A string argument naming a table the workspace has not got is refused", {tag: ["@realizes:powerpack.dialogs.add-new-column", "@realizes:GROK-20931", "@known-failure:27", "@realizes:public#4404"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens a table \"funcs\" with:", () => openTableOf(page, "funcs", [["x"],["1"],["2"],["3"]]), [["x"],["1"],["2"],["3"]]);
    await knownFailure(async () => {
      await session.step(25, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
      await session.step(26, "And user types \"RowCount(\\\"nope\\\")\" into formula editor", () => typeInto(page, "RowCount(\"nope\")", el("formula editor")));
      await session.step(27, "Then formula error should contain text \"nope\"", () => shouldContainText(page, el("formula error"), "nope"));
      await session.step(28, "And OK button in \"Add New Column\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Add New Column\" dialog"), "disabled"));
      await session.step(29, "And no errors should have been logged", () => noErrors(page));
    }, 27);
  });
  test("On auto, the Type choice of Avg names the type of the column OK adds", {tag: ["@realizes:powerpack.dialogs.add-new-column", "@realizes:GROK-20931", "@known-failure:35", "@realizes:public#4406"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens a table \"funcs\" with:", () => openTableOf(page, "funcs", [["x"],["1"],["2"],["3"]]), [["x"],["1"],["2"],["3"]]);
    await knownFailure(async () => {
      await session.step(33, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
      await session.step(34, "And user types \"Avg($[x])\" into formula editor", () => typeInto(page, "Avg($[x])", el("formula editor")));
      await session.step(35, "Then column type input should have value \"auto (double)\"", () => shouldHaveValue(page, el("column type input"), "auto (double)"));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
    }, 35);
  });
});
