/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/pmx-studio-teaser.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {clickOpensStage, clickSettles, enterSettles, guideBudget, openOverview, pointAtMenu, presenting, studiesShown} from '../../bindings/pmx.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, dragSliderTo, enterInto, scrollIntoView, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("PMX Studio", () => {
  const session = feature(test, "features/guides/pmx-studio-teaser.feature", import.meta.url);
  test("An efficacy study from its design to its analysis", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given the guide may take up to 8 minutes", () => guideBudget(page, 8));
    await session.step(9, "And user is logged in", () => loggedIn(page));
    await session.step(10, "And simple mode is off", () => simpleModeOff(page));
    await session.step(11, "And the shell is in presentation mode without tooltips", () => presenting(page));
    await session.step(12, "And user opens PMX Studio on its Overview", () => openOverview(page));
    await session.step(14, "Then the lab's studies should be on the Overview", () => studiesShown(page), undefined, "Every tumor efficacy study of the lab, from its design to its analysis");
    await session.step(16, "When user clicks on demo plan row and the stage opens", () => clickOpensStage(page, el("demo plan row")), undefined, "Plan: group sizes, duration and power, predicted from historical control studies");
    await session.step(18, "And user enters \"1200\" into \"Endpoint volume\" input and PMX Studio settles", () => enterSettles(page, "1200", el("\"Endpoint volume\" input")), undefined, "Every assumption is live: lower the endpoint volume and the design follows");
    await session.step(20, "When user clicks on monitor step and the stage opens", () => clickOpensStage(page, el("monitor step")), undefined, "Monitor: the running study against what its plan predicted");
    await session.step(22, "And user clicks on \"Follow an example study\" button and PMX Studio settles", () => clickSettles(page, el("\"Follow an example study\" button")), undefined, "Monitor: the running study against what its plan predicted");
    await session.step(24, "Then \"Needs attention\" section should be visible", () => shouldBe(page, el("\"Needs attention\" section"), "visible"), undefined, "Halfway through: each group against its forecast, and when the study should end");
    await session.step(26, "When user clicks on \"Measurements\" tab", () => clickOn(page, el("\"Measurements\" tab")), undefined, "For every animal, the range its next measurement should fall in");
    await session.step(28, "And user picks \"Flag this animal…\" from the context menu of first sheet row", () => pickFromContextMenu(page, "Flag this animal…", el("first sheet row")), undefined, "Any animal can be flagged for review, with its reason");
    await session.step(30, "And user enters \"Ulcerated tumor: check at the next session\" into \"Reason\" input", () => enterInto(page, "Ulcerated tumor: check at the next session", el("\"Reason\" input")), undefined, "Any animal can be flagged for review, with its reason");
    await session.step(32, "And user clicks on \"OK\" button and PMX Studio settles", () => clickSettles(page, el("\"OK\" button")), undefined, "Flags and exclusions are logged with who, when and why, and travel with the study");
    await session.step(34, "When user clicks on analyze header step and the stage opens", () => clickOpensStage(page, el("analyze header step")), undefined, "Analyze: the finished study, every arm against the control");
    await session.step(36, "And user clicks on \"Analyze an example study\" button and PMX Studio settles", () => clickSettles(page, el("\"Analyze an example study\" button")), undefined, "Analyze: the finished study, every arm against the control");
    await session.step(38, "When user clicks on each animal tab and PMX Studio settles", () => clickSettles(page, el("each animal tab")), undefined, "Each animal against its digital twin: the same tumor, predicted untreated");
    await session.step(40, "When user clicks on results tab", () => clickOn(page, el("results tab")), undefined, "Time to the endpoint, as measured");
    await session.step(42, "And user scrolls km chart into view", () => scrollIntoView(page, el("km chart")), undefined, "Time to the endpoint, as measured");
    await session.step(44, "And user drags the slider of km slider to 40", () => dragSliderTo(page, el("km slider"), 40), undefined, "… and projected past the last measurement, marked as a projection");
    await session.step(46, "When user points at \"Analyze > Combination synergy…\" in the top menu", () => pointAtMenu(page, "Analyze > Combination synergy…"), undefined, "Coming next: combination therapies, A + B against A and B, synergy or antagonism");
    await session.step(48, "When user clicks on overview step and the stage opens", () => clickOpensStage(page, el("overview step")), undefined, "The whole portfolio in one place: planned, running and finished");
  });
});
