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
import {clickOpensStage, clickSettles, defaultLayouts, enterSettles, noTooltips, scrollTo, studiesListed} from '../../bindings/pmx.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, expand, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Explore PMX Studio", () => {
  const session = feature(test, "features/guides/explore-pmx-studio.feature", import.meta.url);
  test.skip("Plan a study, follow it and analyze it", {tag: ["@skip", "@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And simple mode is off", () => simpleModeOff(page));
    await session.step(13, "And the stages show their default layouts", () => defaultLayouts(page));
    await session.step(14, "And the shell shows no tooltips", () => noTooltips(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(17, "When user expands Apps tree node inside browse tree", () => expand(page, el("Apps tree node inside browse tree")), undefined, "Expand Apps in the Browse panel");
    await session.step(19, "And user clicks on Apps---PMX-Studio tree node inside browse tree", () => clickOn(page, el("Apps---PMX-Studio tree node inside browse tree")), undefined, "Open PMX Studio: the Overview lists your studies");
    await session.step(20, "Then the studies should be listed", () => studiesListed(page));
    await session.step(22, "When user clicks on demo plan row and the stage opens", () => clickOpensStage(page, el("demo plan row")), undefined, "Open a planned study");
    await session.step(24, "Then \"Predicted tumor growth\" section should be visible", () => shouldBe(page, el("\"Predicted tumor growth\" section"), "visible"), undefined, "The plan: each group's predicted growth, from historical control studies");
    await session.step(26, "When user enters \"1200\" into \"Endpoint volume\" input and PMX Studio settles", () => enterSettles(page, "1200", el("\"Endpoint volume\" input")), undefined, "Change the endpoint volume: the design follows");
    await session.step(28, "And user scrolls the page to study calendar", () => scrollTo(page, el("study calendar")), undefined, "Scroll to the study calendar");
    await session.step(30, "When user clicks on monitor step and the stage opens", () => clickOpensStage(page, el("monitor step")), undefined, "Go to Monitor from the header");
    await session.step(32, "And user clicks on \"Follow an example study\" button and PMX Studio settles", () => clickSettles(page, el("\"Follow an example study\" button")), undefined, "Follow an example study");
    await session.step(34, "Then \"Needs attention\" section should be visible", () => shouldBe(page, el("\"Needs attention\" section"), "visible"), undefined, "The study today: each group against its forecast, and who needs attention");
    await session.step(36, "When user clicks on \"Measurements\" tab", () => clickOn(page, el("\"Measurements\" tab")), undefined, "See each animal's measurements and the range expected next");
    await session.step(38, "When user clicks on \"Welfare\" tab and PMX Studio settles", () => clickSettles(page, el("\"Welfare\" tab")), undefined, "Check body weight and condition");
    await session.step(40, "When user clicks on analyze header step and the stage opens", () => clickOpensStage(page, el("analyze header step")), undefined, "Go to Analyze from the header: the study you follow, analyzed as of today");
    await session.step(42, "Then \"Experimental groups\" section should be visible", () => shouldBe(page, el("\"Experimental groups\" section"), "visible"), undefined, "Each treated group against the control");
    await session.step(44, "When user clicks on each animal tab and PMX Studio settles", () => clickSettles(page, el("each animal tab")), undefined, "Compare each animal with its untreated twin");
    await session.step(46, "Then \"Each animal against its untreated twin\" section should be visible", () => shouldBe(page, el("\"Each animal against its untreated twin\" section"), "visible"), undefined, "The gap between each animal and its twin is the effect of the treatment");
    await session.step(48, "When user clicks on \"Tolerability\" tab and PMX Studio settles", () => clickSettles(page, el("\"Tolerability\" tab")), undefined, "Check tolerability");
    await session.step(50, "When user clicks on \"Against the plan\" tab and PMX Studio settles", () => clickSettles(page, el("\"Against the plan\" tab")), undefined, "Compare the study with its plan");
    await session.step(52, "Then \"Plan vs observed\" section should be visible", () => shouldBe(page, el("\"Plan vs observed\" section"), "visible"), undefined, "What the study found against what its plan predicted");
    await session.step(54, "When user clicks on overview step and the stage opens", () => clickOpensStage(page, el("overview step")), undefined, "Back to the Overview");
  });
});
