/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demo/platform/filters.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.filter-builder, u2.filter-query-input, u2.tabs]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {openDemoPage} from '../../../bindings/demo.js';
import {tallPx} from '@datagrok-libraries/bdd/bindings/common/geometry';
import {clickOn, enterInto, selectIn, shouldHaveText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {shouldShowTab} from '@datagrok-libraries/bdd/bindings/common/u2-steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Filter builder", () => {
  const session = feature(test, "features/demo/platform/filters.feature", import.meta.url);
  test("A condition built in the builder filters the dataframe on Apply", {tag: ["@demo", "@realizes:u2.filter-builder", "@realizes:u2.filter-query-input", "@realizes:u2.tabs"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Filters\" demo page", () => openDemoPage(page, "Filters"));
    enter(page, "U2 Demo");
    await session.step(11, "Then \"Entity type\" tab should be 34 px tall", () => tallPx(page, el("\"Entity type\" tab"), 34));
    await session.step(12, "When user clicks on \"Add condition\" button", () => clickOn(page, el("\"Add condition\" button")));
    await session.step(13, "And user selects \"city\" in property of filter row", () => selectIn(page, "city", el("property of filter row")));
    await session.step(14, "And user selects \"equals\" in operator of filter row", () => selectIn(page, "equals", el("operator of filter row")));
    await session.step(15, "And user enters \"Basel\" into value of filter row", () => enterInto(page, "Basel", el("value of filter row")));
    await session.step(16, "Then Query input should have value 'city = \"Basel\"'", () => shouldHaveValue(page, el("Query input"), "city = \"Basel\""));
    await session.step(17, "When user clicks on Apply button", () => clickOn(page, el("Apply button")));
    await session.step(18, "Then value of rows readout should have text \"2 of 6\"", () => shouldHaveText(page, el("value of rows readout"), "2 of 6"));
    await session.step(19, "When user clicks on \"Entity type\" tab", () => clickOn(page, el("\"Entity type\" tab")));
    await session.step(20, "Then tabs should show the \"Entity type\" tab", () => shouldShowTab(page, el("tabs"), "Entity type"));
    await session.step(21, "When user clicks on DataFrame tab", () => clickOn(page, el("DataFrame tab")));
    await session.step(22, "Then Query input should have value 'city = \"Basel\"'", () => shouldHaveValue(page, el("Query input"), "city = \"Basel\""));
    await session.step(23, "And value of rows readout should have text \"2 of 6\"", () => shouldHaveText(page, el("value of rows readout"), "2 of 6"));
  });
  test("A condition on the users table counts its rows", {tag: ["@demo", "@realizes:u2.filter-builder", "@realizes:u2.filter-query-input", "@realizes:u2.tabs"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Filters\" demo page", () => openDemoPage(page, "Filters"));
    enter(page, "U2 Demo");
    await session.step(26, "When user clicks on \"Domain table\" tab", () => clickOn(page, el("\"Domain table\" tab")));
    await session.step(27, "And user clicks on Open button", () => clickOn(page, el("Open button")));
    await session.step(28, "And user clicks on \"Add condition\" button", () => clickOn(page, el("\"Add condition\" button")));
    await session.step(29, "And user selects \"Login\" in property of filter row", () => selectIn(page, "Login", el("property of filter row")));
    await session.step(30, "And user selects \"equals\" in operator of filter row", () => selectIn(page, "equals", el("operator of filter row")));
    await session.step(31, "And user enters \"admin\" into value of filter row", () => enterInto(page, "admin", el("value of filter row")));
    await session.step(32, "Then Query input should have value 'login = \"admin\"'", () => shouldHaveValue(page, el("Query input"), "login = \"admin\""));
    await session.step(33, "When user clicks on Count button", () => clickOn(page, el("Count button")));
    await session.step(34, "Then value of rows readout should have text \"1\"", () => shouldHaveText(page, el("value of rows readout"), "1"));
  });
});
