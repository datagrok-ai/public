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
import {finishedRunning, guideBudget, pointAtMenu, presenting, stageOpen, studiesListed} from '../../bindings/pmx.js';
import {openApp} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, dragSliderTo, enterInto, scrollIntoView, shouldContainText, waitSeconds} from '@datagrok-libraries/bdd/bindings/common/steps';
import {simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("PMX Studio in two minutes", () => {
  const session = feature(test, "features/guides/pmx-studio-teaser.feature", import.meta.url);
  test("One study from its plan to its analysis", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given the guide may take up to 8 minutes", () => guideBudget(page, 8));
    await session.step(9, "And user is logged in", () => loggedIn(page));
    await session.step(10, "And simple mode is off", () => simpleModeOff(page));
    await session.step(11, "And the shell is in presentation mode without tooltips", () => presenting(page));
    await session.step(13, "And user opens PMX Studio", () => openApp(page), undefined, "PMX Studio: the lab's studies, planned, running and finished");
    await session.step(14, "Then the studies should be listed", () => studiesListed(page));
    await session.step(16, "And user waits 3 seconds", () => waitSeconds(page, 3), undefined, "PMX Studio: the lab's studies, planned, running and finished");
    await session.step(18, "When user clicks on demo plan row", () => clickOn(page, el("demo plan row")), undefined, "The plan of an MC38 study, sized from the atlas of historical controls");
    await session.step(19, "Then the plan stage should be open", () => stageOpen(page, "plan"));
    await session.step(21, "And user waits 2 seconds", () => waitSeconds(page, 2), undefined, "The plan of an MC38 study, sized from the atlas of historical controls");
    await session.step(23, "When user enters \"1200\" into \"Endpoint volume\" input", () => enterInto(page, "1200", el("\"Endpoint volume\" input")), undefined, "Change an assumption: the design follows at once");
    await session.step(24, "Then the study should have finished running", () => finishedRunning(page));
    await session.step(26, "And user waits 2 seconds", () => waitSeconds(page, 2), undefined, "Change an assumption: the design follows at once");
    await session.step(28, "When user clicks on overview step", () => clickOn(page, el("overview step")), undefined, "The same study, halfway through");
    await session.step(30, "And user clicks on demo followed row", () => clickOn(page, el("demo followed row")), undefined, "The same study, halfway through");
    await session.step(31, "Then the monitor stage should be open", () => stageOpen(page, "monitor"));
    await session.step(33, "And user waits 3 seconds", () => waitSeconds(page, 3), undefined, "Monitor: who needs attention today, and when the study should end");
    await session.step(35, "When user picks \"Flag this animal…\" from the context menu of first attention row", () => pickFromContextMenu(page, "Flag this animal…", el("first attention row")), undefined, "Flag an animal from its row");
    await session.step(37, "And user enters \"Ulcerated tumor: check at the next session\" into \"Reason\" input", () => enterInto(page, "Ulcerated tumor: check at the next session", el("\"Reason\" input")), undefined, "With its reason: it goes in the study's log");
    await session.step(39, "And user clicks on \"OK\" button", () => clickOn(page, el("\"OK\" button")), undefined, "With its reason: it goes in the study's log");
    await session.step(40, "Then the study should have finished running", () => finishedRunning(page));
    await session.step(42, "And user waits 2 seconds", () => waitSeconds(page, 2), undefined, "The forecasts are recomputed with it");
    await session.step(44, "When user clicks on overview step", () => clickOn(page, el("overview step")), undefined, "The study with all its data");
    await session.step(46, "And user clicks on demo finished row", () => clickOn(page, el("demo finished row")), undefined, "The study with all its data");
    await session.step(47, "Then the analyze stage should be open", () => stageOpen(page, "analyze"));
    await session.step(49, "And user waits 2 seconds", () => waitSeconds(page, 2), undefined, "Analyze: every arm against the control");
    await session.step(51, "When user clicks on each animal tab", () => clickOn(page, el("each animal tab")), undefined, "Each animal against its digital twin, untreated");
    await session.step(52, "Then the study should have finished running", () => finishedRunning(page));
    await session.step(54, "And user waits 3 seconds", () => waitSeconds(page, 3), undefined, "Each animal against its digital twin, untreated");
    await session.step(56, "When user clicks on results tab", () => clickOn(page, el("results tab")), undefined, "Survival to the endpoint");
    await session.step(58, "And user scrolls km chart into view", () => scrollIntoView(page, el("km chart")), undefined, "Survival to the endpoint");
    await session.step(60, "And user drags the slider of km slider to 40", () => dragSliderTo(page, el("km slider"), 40), undefined, "Follow the study past its last measurement");
    await session.step(62, "Then \"Time to endpoint\" section should contain text \"projected\"", () => shouldContainText(page, el("\"Time to endpoint\" section"), "projected"), undefined, "A projection, marked as such");
    await session.step(64, "And user waits 2 seconds", () => waitSeconds(page, 2), undefined, "A projection, marked as such");
    await session.step(66, "When user points at \"Analyze > Combination synergy…\" in the top menu", () => pointAtMenu(page, "Analyze > Combination synergy…"), undefined, "Next: combination therapies. A, B and A+B: synergy or antagonism");
    await session.step(68, "When user clicks on overview step", () => clickOn(page, el("overview step")), undefined, "Back on the Overview: planned, followed and finished, in one place");
    await session.step(69, "Then the studies should be listed", () => studiesListed(page));
    await session.step(71, "And user waits 3 seconds", () => waitSeconds(page, 3), undefined, "Back on the Overview: planned, followed and finished, in one place");
  });
});
