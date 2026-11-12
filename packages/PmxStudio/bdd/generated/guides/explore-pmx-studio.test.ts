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
import {clickOpensStage, clickSettles, defaultLayouts, enterSettles, guideBudget, noSavedPlan, noTooltips, planIsSaved, scrollTo} from '../../bindings/pmx.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, expand, shouldBe, uploadThrough} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Explore PMX Studio", () => {
  const session = feature(test, "features/guides/explore-pmx-studio.feature", import.meta.url);
  test("Plan a study, follow it and analyze it", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given the guide may take up to 15 minutes", () => guideBudget(page, 15));
    await session.step(10, "And user is logged in", () => loggedIn(page));
    await session.step(11, "And simple mode is off", () => simpleModeOff(page));
    await session.step(12, "And the stages show their default layouts", () => defaultLayouts(page));
    await session.step(13, "And no plan named \"MC38 guide plan\" is saved", () => noSavedPlan(page, "MC38 guide plan"));
    await session.step(14, "And the shell shows no tooltips", () => noTooltips(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(17, "When user expands Apps tree node inside browse tree", () => expand(page, el("Apps tree node inside browse tree")), undefined, "Expand Apps in the Browse panel");
    await session.step(19, "And user clicks on Apps---PMX-Studio tree node inside browse tree", () => clickOn(page, el("Apps---PMX-Studio tree node inside browse tree")), undefined, "Open PMX Studio: the Overview lists your studies");
    await session.step(20, "Then open a study strip should be visible", () => shouldBe(page, el("open a study strip"), "visible"));
    await session.step(22, "When user clicks on plan card and the stage opens", () => clickOpensStage(page, el("plan card")), undefined, "Start a new study in Plan: the tumor model comes first");
    await session.step(24, "And user clicks on \"Add treatment arm\" button and PMX Studio settles", () => clickSettles(page, el("\"Add treatment arm\" button")), undefined, "Add the treatment arms the study compares with its control");
    await session.step(25, "And user clicks on \"Add treatment arm\" button and PMX Studio settles", () => clickSettles(page, el("\"Add treatment arm\" button")));
    await session.step(26, "And user clicks on \"Add treatment arm\" button and PMX Studio settles", () => clickSettles(page, el("\"Add treatment arm\" button")));
    await session.step(28, "Then \"Predicted tumor growth\" section should be visible", () => shouldBe(page, el("\"Predicted tumor growth\" section"), "visible"), undefined, "The plan: each group's predicted growth, from historical control studies");
    await session.step(30, "When user enters \"1200\" into \"Endpoint volume\" input and PMX Studio settles", () => enterSettles(page, "1200", el("\"Endpoint volume\" input")), undefined, "Change the endpoint volume: the design follows");
    await session.step(32, "And user scrolls the page to study calendar", () => scrollTo(page, el("study calendar")), undefined, "Scroll to the study calendar");
    await session.step(34, "When user clicks on \"Save plan\" button", () => clickOn(page, el("\"Save plan\" button")), undefined, "Save the plan");
    await session.step(35, "And user enters \"MC38 guide plan\" into \"Name\" input", () => enterInto(page, "MC38 guide plan", el("\"Name\" input")));
    await session.step(36, "And user clicks on \"Save\" button in \"Save plan\" dialog", () => clickOn(page, el("\"Save\" button in \"Save plan\" dialog")));
    await session.step(37, "Then a plan named \"MC38 guide plan\" should be saved", () => planIsSaved(page, "MC38 guide plan"));
    await session.step(39, "When user clicks on overview step and the stage opens", () => clickOpensStage(page, el("overview step")), undefined, "Back on the Overview, open a study: its plan, then its measurements so far");
    await session.step(40, "And user clicks on open a study strip", () => clickOn(page, el("open a study strip")));
    await session.step(41, "And user clicks on saved plans link", () => clickOn(page, el("saved plans link")));
    await session.step(42, "And user clicks on guide plan menu item", () => clickOn(page, el("guide plan menu item")));
    await session.step(44, "And user uploads \"../files/example/mc38-ld4172-anti-pd1.csv\" through second file zone", () => uploadThrough(page, "../files/example/mc38-ld4172-anti-pd1.csv", el("second file zone")), undefined, "The measurements: here an example table, a published MC38 study");
    await session.step(46, "And user clicks on still running choice", () => clickOn(page, el("still running choice")), undefined, "Still running: follow it in Monitor");
    await session.step(47, "And user clicks on \"Open in Monitor\" button in \"Open a study\" dialog and the stage opens", () => clickOpensStage(page, el("\"Open in Monitor\" button in \"Open a study\" dialog")));
    await session.step(49, "Then \"Needs attention\" section should be visible", () => shouldBe(page, el("\"Needs attention\" section"), "visible"), undefined, "The study today: each group against its forecast, and who needs attention");
    await session.step(51, "When user clicks on \"Measurements\" tab", () => clickOn(page, el("\"Measurements\" tab")), undefined, "See each animal's measurements and the range expected next");
    await session.step(53, "When user clicks on \"Welfare\" tab and PMX Studio settles", () => clickSettles(page, el("\"Welfare\" tab")), undefined, "Check body weight and condition");
    await session.step(55, "When user clicks on analyze header step and the stage opens", () => clickOpensStage(page, el("analyze header step")), undefined, "Go to Analyze from the header: the study you follow, analyzed as of today");
    await session.step(57, "Then \"Experimental groups\" section should be visible", () => shouldBe(page, el("\"Experimental groups\" section"), "visible"), undefined, "Each treated group against the control");
    await session.step(59, "When user clicks on each animal tab and PMX Studio settles", () => clickSettles(page, el("each animal tab")), undefined, "Compare each animal with its untreated twin");
    await session.step(61, "Then \"Each animal against its untreated twin\" section should be visible", () => shouldBe(page, el("\"Each animal against its untreated twin\" section"), "visible"), undefined, "The gap between each animal and its twin is the effect of the treatment");
    await session.step(63, "When user clicks on \"Tolerability\" tab and PMX Studio settles", () => clickSettles(page, el("\"Tolerability\" tab")), undefined, "Check tolerability");
    await session.step(65, "When user clicks on \"Against the plan\" tab and PMX Studio settles", () => clickSettles(page, el("\"Against the plan\" tab")), undefined, "Compare the study with its plan");
    await session.step(67, "Then \"Plan vs observed\" section should be visible", () => shouldBe(page, el("\"Plan vs observed\" section"), "visible"), undefined, "What the study found against what its plan predicted");
    await session.step(69, "When user clicks on overview step and the stage opens", () => clickOpensStage(page, el("overview step")), undefined, "Back to the Overview");
  });
});
