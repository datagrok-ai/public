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
import {clickOpensStage, clickSettles, defaultLayouts, demoOpen, dragSettles, guideBudget, openOverview, pointOut, presenting, scrolledTo, studiesShown} from '../../bindings/pmx.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, dragSliderTo, enterInto, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("PMX Studio", () => {
  const session = feature(test, "features/guides/pmx-studio-teaser.feature", import.meta.url);
  test("An efficacy study from its design to its analysis", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given the guide may take up to 10 minutes", () => guideBudget(page, 10));
    await session.step(9, "And user is logged in", () => loggedIn(page));
    await session.step(10, "And simple mode is off", () => simpleModeOff(page));
    await session.step(11, "And the stages show their default layouts", () => defaultLayouts(page));
    await session.step(12, "And the shell is in presentation mode without tooltips", () => presenting(page));
    await session.step(13, "And user opens PMX Studio on its Overview", () => openOverview(page));
    await session.step(14, "And the demo study is open in Monitor and Analyze", () => demoOpen(page));
    await session.step(16, "Then the lab's studies should be on the Overview", () => studiesShown(page), undefined, "Every tumor efficacy study of the lab, from its design to its analysis");
    await session.step(18, "When user points out your studies", () => pointOut(page, el("your studies")), undefined, "Planned, running and finished studies side by side, each one a click away");
    await session.step(20, "When user clicks on demo plan row and the stage opens", () => clickOpensStage(page, el("demo plan row")), undefined, "Plan: group sizes, duration and power, predicted from historical control studies");
    await session.step(22, "And user points out tumor model section", () => pointOut(page, el("tumor model section")), undefined, "The design starts from the tumor model: its growth and variability come from past control groups");
    await session.step(24, "And user drags the slider of \"Endpoint volume\" input to 1200 and PMX Studio settles", () => dragSettles(page, el("\"Endpoint volume\" input"), 1200), undefined, "Every assumption is live: lower the endpoint volume and the design follows");
    await session.step(26, "And user points out plan figures", () => pointOut(page, el("plan figures")), undefined, "Animals per group, study length and power, recomputed at once");
    await session.step(28, "And user points out treatment arms section", () => pointOut(page, el("treatment arms section")), undefined, "The arms to compare, each with its dose and schedule");
    await session.step(30, "When user clicks on monitor step and the stage opens", () => clickOpensStage(page, el("monitor step")), undefined, "Monitor: the running study against what its plan predicted");
    await session.step(32, "Then \"Needs attention\" section should be visible", () => shouldBe(page, el("\"Needs attention\" section"), "visible"), undefined, "Each group against its forecast, who needs attention today, and when the study should end");
    await session.step(34, "When user clicks on \"Measurements\" tab", () => clickOn(page, el("\"Measurements\" tab")), undefined, "For every animal, the range its next measurement should fall in");
    await session.step(36, "And user picks \"Flag this animal…\" from the context menu of first sheet row", () => pickFromContextMenu(page, "Flag this animal…", el("first sheet row")), undefined, "Any animal can be flagged for review, with its reason");
    await session.step(38, "And user enters \"Ulcerated tumor: check at the next session\" into \"Reason\" input", () => enterInto(page, "Ulcerated tumor: check at the next session", el("\"Reason\" input")), undefined, "Any animal can be flagged for review, with its reason");
    await session.step(40, "And user clicks on \"OK\" button and PMX Studio settles", () => clickSettles(page, el("\"OK\" button")), undefined, "Any animal can be flagged for review, with its reason");
    await session.step(42, "And user points out first sheet row", () => pointOut(page, el("first sheet row")), undefined, "Flagged: the animal carries its mark, and its rows are coloured in the study table");
    await session.step(44, "When user clicks on \"Welfare\" tab and PMX Studio settles", () => clickSettles(page, el("\"Welfare\" tab")), undefined, "Welfare: body weight and condition of every animal, against the humane endpoints");
    await session.step(46, "When user clicks on \"Decisions\" tab and PMX Studio settles", () => clickSettles(page, el("\"Decisions\" tab")), undefined, "Decisions: flags and exclusions logged with who, when and why; they travel with the study");
    await session.step(48, "When user clicks on analyze header step and the stage opens", () => clickOpensStage(page, el("analyze header step")), undefined, "Analyze: the finished study, every arm against the control");
    await session.step(50, "And user points out experimental groups section", () => pointOut(page, el("experimental groups section")), undefined, "Each group's growth and its effect against the control, with its uncertainty");
    await session.step(51, "And the page is scrolled to km chart", () => scrolledTo(page, el("km chart")));
    await session.step(53, "And user drags the slider of km slider to 40", () => dragSliderTo(page, el("km slider"), 40), undefined, "Time to the endpoint as measured, then projected past the last measurement, marked as a projection");
    await session.step(54, "And the page is scrolled to each animal tab", () => scrolledTo(page, el("each animal tab")));
    await session.step(56, "When user clicks on each animal tab and PMX Studio settles", () => clickSettles(page, el("each animal tab")), undefined, "Each animal: its response, from the best to the worst, and its tumor over time");
    await session.step(58, "And user points out waterfall chart", () => pointOut(page, el("waterfall chart")), undefined, "Each animal: its response, from the best to the worst, and its tumor over time");
    await session.step(59, "And the page is scrolled to twin section", () => scrolledTo(page, el("twin section")));
    await session.step(61, "And user clicks on third animal row and PMX Studio settles", () => clickSettles(page, el("third animal row")), undefined, "Its digital twin: the same tumor, predicted as if untreated");
    await session.step(63, "And user points out twin effect line", () => pointOut(page, el("twin effect line")), undefined, "The gap between the animal and its twin is the effect of the treatment, animal by animal");
    await session.step(64, "And the page is scrolled to tolerability tab", () => scrolledTo(page, el("tolerability tab")));
    await session.step(66, "When user clicks on \"Tolerability\" tab and PMX Studio settles", () => clickSettles(page, el("\"Tolerability\" tab")), undefined, "Tolerability: body weight loss and deaths in every arm");
    await session.step(68, "When user clicks on \"Against the plan\" tab and PMX Studio settles", () => clickSettles(page, el("\"Against the plan\" tab")), undefined, "Against the plan: what the study found against what its design predicted");
    await session.step(70, "When user clicks on overview step and the stage opens", () => clickOpensStage(page, el("overview step")), undefined, "The whole portfolio in one place: planned, running and finished");
  });
});
