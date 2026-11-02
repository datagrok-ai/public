/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/explore-pmx-studio.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {finishedRunning} from '../../bindings/pmx.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, expand, selectIn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Explore PMX Studio", () => {
  const session = feature(test, "features/guides/explore-pmx-studio.feature", import.meta.url);
  test("Plan a study and analyze an example", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(7, "Given user is logged in", () => loggedIn(page));
    await session.step(8, "And simple mode is off", () => simpleModeOff(page));
    await session.step(9, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(11, "When user expands Apps tree node inside browse tree", () => expand(page, el("Apps tree node inside browse tree")), undefined, "Expand Apps in the Browse panel");
    await session.step(13, "And user clicks on Apps---PMX-Studio tree node inside browse tree", () => clickOn(page, el("Apps---PMX-Studio tree node inside browse tree")), undefined, "Open PMX Studio");
    await session.step(14, "Then the study should have finished running", () => finishedRunning(page));
    await session.step(15, "And \"Experimental groups\" section should contain text \"Control\"", () => shouldContainText(page, el("\"Experimental groups\" section"), "Control"));
    await session.step(17, "When user selects \"CT26\" in \"Cell line\" input", () => selectIn(page, "CT26", el("\"Cell line\" input")), undefined, "Choose the CT26 cell line");
    await session.step(19, "And user clicks on form run button", () => clickOn(page, el("form run button")), undefined, "Run the plan for the new cell line");
    await session.step(20, "Then the study should have finished running", () => finishedRunning(page));
    await session.step(22, "When user clicks on \"Calendar\" tab", () => clickOn(page, el("\"Calendar\" tab")), undefined, "See the study calendar");
    await session.step(23, "Then \"Study calendar\" section should be visible", () => shouldBe(page, el("\"Study calendar\" section"), "visible"));
    await session.step(25, "When user clicks on \"Trade-offs\" tab", () => clickOn(page, el("\"Trade-offs\" tab")), undefined, "See how power changes with the design");
    await session.step(26, "Then \"Power by group size\" section should be visible", () => shouldBe(page, el("\"Power by group size\" section"), "visible"));
    await session.step(28, "When user clicks on \"Protocol\" tab", () => clickOn(page, el("\"Protocol\" tab")), undefined, "Read the protocol, ready to paste");
    await session.step(29, "Then \"Study protocol\" section should contain text \"CT26\"", () => shouldContainText(page, el("\"Study protocol\" section"), "CT26"));
    await session.step(31, "When user clicks on \"3 Analyze\" button", () => clickOn(page, el("\"3 Analyze\" button")), undefined, "Go to the Analyze stage");
    await session.step(33, "And user clicks on \"Analyze an example study\" button", () => clickOn(page, el("\"Analyze an example study\" button")), undefined, "Analyze an example study");
    await session.step(34, "Then the study should have finished running", () => finishedRunning(page));
    await session.step(35, "And \"Observed tumor growth\" section should be visible", () => shouldBe(page, el("\"Observed tumor growth\" section"), "visible"));
    await session.step(37, "When user clicks on \"Each animal\" tab", () => clickOn(page, el("\"Each animal\" tab")), undefined, "Compare each animal with its untreated twin");
    await session.step(38, "Then \"Each animal against its untreated twin\" section should be visible", () => shouldBe(page, el("\"Each animal against its untreated twin\" section"), "visible"));
    await session.step(40, "When user clicks on \"Tolerability\" tab", () => clickOn(page, el("\"Tolerability\" tab")), undefined, "Check tolerability");
    await session.step(41, "Then \"Weight loss\" section should be visible", () => shouldBe(page, el("\"Weight loss\" section"), "visible"));
    await session.step(43, "When user clicks on \"Plan vs observed\" tab", () => clickOn(page, el("\"Plan vs observed\" tab")), undefined, "Compare the study with its plan");
    await session.step(44, "Then \"Effect of each treatment\" section should be visible", () => shouldBe(page, el("\"Effect of each treatment\" section"), "visible"));
    await session.step(46, "When user clicks on \"Statistics\" tab", () => clickOn(page, el("\"Statistics\" tab")), undefined, "See the statistics");
    await session.step(47, "Then \"Results summary\" section should be visible", () => shouldBe(page, el("\"Results summary\" section"), "visible"));
  });
});
