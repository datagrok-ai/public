/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/filter-builder/filter-builder.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
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
import {clickOn, enterInto, selectIn, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterIsExactly, filterPasses, filterPassesAll, openFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickPanelMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The filter panel's Filter Builder filters the table by conditions built in rows or typed as a query", () => {
  const session = feature(test, "features/filter-builder/filter-builder.feature", import.meta.url);
  test("The filter panel's Filter Builder filters the table by conditions built in rows or typed as a query", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And user opens a table \"people\" with:", () => openTableOf(page, "people", [["name","v"],["Ann","1"],["Bob","2"],["Cid","3"],["Dan","4"]]), [["name","v"],["Ann","1"],["Bob","2"],["Cid","3"],["Dan","4"]]);
    await session.step(28, "And user opens the filter panel", () => openFilterPanel(page));
    await session.step(29, "And user picks \"Add Filter | Filter Builder\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Filter Builder"));
    await run.scenario("A condition built in a row filters the table, and removing the row lifts the filter", async () => {
      await session.step(32, "Then \"Add condition\" link in filter builder should be visible", () => shouldBe(page, el("\"Add condition\" link in filter builder"), "visible"));
      await session.step(33, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(34, "When user clicks on \"Add condition\" link in filter builder", () => clickOn(page, el("\"Add condition\" link in filter builder")));
      await session.step(35, "And user selects \"v\" in property of first filter row", () => selectIn(page, "v", el("property of first filter row")));
      await session.step(36, "And user selects \"greater than\" in operator of first filter row", () => selectIn(page, "greater than", el("operator of first filter row")));
      await session.step(37, "And user enters \"2\" into value of first filter row", () => enterInto(page, "2", el("value of first filter row")));
      await session.step(38, "Then 2 rows should pass the filter", () => filterPasses(page, 2));
      await session.step(39, "And the filter should pass exactly the rows where \"v\" is between 3 and 4", () => filterIsExactly(page, "v", 3, 4));
      await session.step(40, "When user clicks on remove of first filter row", () => clickOn(page, el("remove of first filter row")));
      await session.step(41, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A query typed into the card's query text filters the table, and the rows show it", async () => {
      await session.step(45, "When user clicks on \"Switch between the builder and the query text\" icon in filter panel", () => clickOn(page, el("\"Switch between the builder and the query text\" icon in filter panel")));
      await session.step(46, "And user enters \"v < 3\" into filter query in filter panel", () => enterInto(page, "v < 3", el("filter query in filter panel")));
      await session.step(47, "Then 2 rows should pass the filter", () => filterPasses(page, 2));
      await session.step(48, "And the filter should pass exactly the rows where \"v\" is between 1 and 2", () => filterIsExactly(page, "v", 1, 2));
      await session.step(49, "When user clicks on \"Switch between the builder and the query text\" icon in filter panel", () => clickOn(page, el("\"Switch between the builder and the query text\" icon in filter panel")));
      await session.step(50, "Then operator of first filter row should have value \"less than\"", () => shouldHaveValue(page, el("operator of first filter row"), "less than"));
      await session.step(51, "And property of first filter row should have value \"v\"", () => shouldHaveValue(page, el("property of first filter row"), "v"));
      await session.step(52, "When user clicks on remove of first filter row", () => clickOn(page, el("remove of first filter row")));
      await session.step(53, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
