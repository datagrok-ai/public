/* ---
generated: features/demo/display/tables.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.table, u2.grid]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openDemoPage} from '../../../bindings/demo.js';
import {clickOn, shouldBe, shouldContainText, shouldHaveRows, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Tables and grids", () => {
  const session = feature(test);
  test("A basic table with selectable rows", {tag: ["@demo", "@realizes:u2.table", "@realizes:u2.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Tables\" demo page", () => openDemoPage(page, "Tables"));
    enter(page, "U2 Demo");
    await test.step("Then table should have 5 rows", () => shouldHaveRows(page, el("table"), 5));
    await test.step("And Caffeine table row should contain text \"28\"", () => shouldContainText(page, el("Caffeine table row"), "28"));
    await test.step("When user clicks on Caffeine table row", () => clickOn(page, el("Caffeine table row")));
    await test.step("Then Caffeine table row should be selected", () => shouldBe(page, el("Caffeine table row"), "selected"));
    await test.step("And value of selectedIndex readout should have text \"1\"", () => shouldHaveText(page, el("value of selectedIndex readout"), "1"));
  });
  test("A virtual grid of 20,000 cells", {tag: ["@demo", "@realizes:u2.table", "@realizes:u2.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Tables\" demo page", () => openDemoPage(page, "Tables"));
    enter(page, "U2 Demo");
    await test.step("When user clicks on \"#7\" item in virtual grid", () => clickOn(page, el("\"#7\" item in virtual grid")));
    await test.step("Then \"#7\" item in virtual grid should be selected", () => shouldBe(page, el("\"#7\" item in virtual grid"), "selected"));
    await test.step("And value of \"selected cell\" readout should have text \"7\"", () => shouldHaveText(page, el("value of \"selected cell\" readout"), "7"));
  });
});
