/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demo/display/app-header.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.app-header]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {openDemoPage} from '../../../bindings/demo.js';
import {clickOn, shouldContainText, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {shouldHaveItemsInTotal} from '@datagrok-libraries/bdd/bindings/common/u2-steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("App header", () => {
  const session = feature(test, "features/demo/display/app-header.feature", import.meta.url);
  test("The bullets are the items, the link items act", {tag: ["@demo", "@realizes:u2.app-header"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(7, "Given user opens the \"App header\" demo page", () => openDemoPage(page, "App header"));
    enter(page, "U2 Demo");
    await session.step(10, "Then \"hit-design\" app header should have 3 items in total", () => shouldHaveItemsInTotal(page, el("\"hit-design\" app header"), 3));
    await session.step(11, "And \"hit-design\" app header should contain text \"Sketch molecules in the molecular spreadsheet\"", () => shouldContainText(page, el("\"hit-design\" app header"), "Sketch molecules in the molecular spreadsheet"));
    await session.step(12, "And MolTrack app header should have 2 items in total", () => shouldHaveItemsInTotal(page, el("MolTrack app header"), 2));
    await session.step(13, "And MolTrack app header should contain text \"Track, analyze, and manage chemical data\"", () => shouldContainText(page, el("MolTrack app header"), "Track, analyze, and manage chemical data"));
    await session.step(14, "When user clicks on \"Open the samples\" link", () => clickOn(page, el("\"Open the samples\" link")));
    await session.step(15, "Then value of actions readout should have text \"1\"", () => shouldHaveText(page, el("value of actions readout"), "1"));
    await session.step(16, "When user clicks on \"Open the help\" link", () => clickOn(page, el("\"Open the help\" link")));
    await session.step(17, "Then value of actions readout should have text \"2\"", () => shouldHaveText(page, el("value of actions readout"), "2"));
  });
});
