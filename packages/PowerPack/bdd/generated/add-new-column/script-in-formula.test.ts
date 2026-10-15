/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/add-new-column/script-in-formula.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [powerpack.dialogs.add-new-column]
--- */
import {test} from '@playwright/test';
import '../../bindings/enrichment.js';
import '../../bindings/home.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {previewComputes} from '../../bindings/add-new-column.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnType, hasColumn, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {openTableOf, scriptOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Add New Column: a script of the user's in a formula", () => {
  const session = feature(test, "features/add-new-column/script-in-formula.feature", import.meta.url);
  test("Add New Column: a script of the user's in a formula", {tag: ["@journey", "@realizes:powerpack.dialogs.add-new-column", "@known-failure:51"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And a script \"BddDouble{time}\" is on the server:", () => scriptOnServer(page, session.text("BddDouble{time}"), "//language: javascript\n//input: double x\n//output: double y\ny = x * 2;"));
    await session.step(23, "And user opens a table \"calc\" with:", () => openTableOf(page, "calc", [["n","f"],["1","1.5"],["2","2.5"],["4","4.5"]]), [["n","f"],["1","1.5"],["2","2.5"],["4","4.5"]]);
    await run.scenario("On a column of fractional numbers the script computes every row", async () => {
      await session.step(30, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
      await session.step(31, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
      await session.step(32, "When user types \"twice\" into column name input", () => typeInto(page, "twice", el("column name input")));
      await session.step(33, "And user types \"BddDouble{time}(${f})\" into formula editor", () => typeInto(page, session.text("BddDouble{time}(${f})"), el("formula editor")));
      await session.step(34, "Then the preview grid should show \"twice\" computed as numbers", () => previewComputes(page, "twice"));
      await session.step(35, "When user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
      await session.step(36, "Then \"Add New Column\" dialog should be hidden", () => shouldBe(page, el("\"Add New Column\" dialog"), "hidden"));
      await session.step(37, "And \"twice\" column should have type \"double\"", () => columnType(page, "twice", "double"));
      await session.step(38, "And the value of \"twice\" column in row 3 should be \"9\"", () => valueInRow(page, "twice", 3, "9"));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("On a column of whole numbers the script computes every row", async () => {
      await session.step(44, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
      await session.step(45, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
      await session.step(46, "When user types \"twice int\" into column name input", () => typeInto(page, "twice int", el("column name input")));
      await session.step(47, "And user types \"BddDouble{time}(${n})\" into formula editor", () => typeInto(page, session.text("BddDouble{time}(${n})"), el("formula editor")));
      await session.step(48, "And user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
      await session.step(49, "Then \"Add New Column\" dialog should be hidden", () => shouldBe(page, el("\"Add New Column\" dialog"), "hidden"));
      await session.step(50, "And the table should have a column \"twice int\"", () => hasColumn(page, "twice int"));
      await session.step(51, "And \"twice int\" column should have type \"double\"", () => columnType(page, "twice int", "double"));
    }, {knownFailure: true, line: 51});
    run.finish();
  });
});
