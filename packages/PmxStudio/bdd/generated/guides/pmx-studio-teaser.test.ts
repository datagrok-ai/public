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
import {clickOpensStage, clickSettles, defaultLayouts, dragSettles, guideBudget, noSavedPlan, openOverview, pointOut, presenting, scrollTo} from '../../bindings/pmx.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, dragSliderTo, enterInto, shouldBe, uploadThrough} from '@datagrok-libraries/bdd/bindings/common/steps';
import {simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("PMX Studio", () => {
  const session = feature(test, "features/guides/pmx-studio-teaser.feature", import.meta.url);
  test("An efficacy study from its design to its analysis", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given the guide may take up to 15 minutes", () => guideBudget(page, 15));
    await session.step(9, "And user is logged in", () => loggedIn(page));
    await session.step(10, "And simple mode is off", () => simpleModeOff(page));
    await session.step(11, "And the stages show their default layouts", () => defaultLayouts(page));
    await session.step(12, "And the shell is in presentation mode without tooltips", () => presenting(page));
    await session.step(13, "And no plan named \"MC38 guide plan\" is saved", () => noSavedPlan(page, "MC38 guide plan"));
    await session.step(14, "And user opens PMX Studio on its Overview", () => openOverview(page));
    await session.step(16, "When user points out your studies", () => pointOut(page, el("your studies")), undefined, "Every efficacy study of the lab in one place: planned, running and finished");
    await session.step(18, "When user clicks on plan card and the stage opens", () => clickOpensStage(page, el("plan card")), undefined, "Plan: the design starts from the tumor model and its historical control studies");
    await session.step(20, "And user points out tumor model section", () => pointOut(page, el("tumor model section")), undefined, "Plan: the design starts from the tumor model and its historical control studies");
    await session.step(22, "And user clicks on \"Add treatment arm\" button and PMX Studio settles", () => clickSettles(page, el("\"Add treatment arm\" button")), undefined, "Three treatment arms, each against the control");
    await session.step(24, "And user clicks on \"Add treatment arm\" button and PMX Studio settles", () => clickSettles(page, el("\"Add treatment arm\" button")), undefined, "Three treatment arms, each against the control");
    await session.step(26, "And user clicks on \"Add treatment arm\" button and PMX Studio settles", () => clickSettles(page, el("\"Add treatment arm\" button")), undefined, "Three treatment arms, each against the control");
    await session.step(28, "And user drags the slider of \"Endpoint volume\" input to 1200 and PMX Studio settles", () => dragSettles(page, el("\"Endpoint volume\" input"), 1200), undefined, "Every assumption is live: lower the endpoint and group sizes, duration and power follow");
    await session.step(30, "And user points out plan figures", () => pointOut(page, el("plan figures")), undefined, "Every assumption is live: lower the endpoint and group sizes, duration and power follow");
    await session.step(32, "And user points out treatment arms section", () => pointOut(page, el("treatment arms section")), undefined, "Every assumption is live: lower the endpoint and group sizes, duration and power follow");
    await session.step(34, "When user clicks on \"Save plan\" button", () => clickOn(page, el("\"Save plan\" button")), undefined, "The plan saved: the study will be followed against it");
    await session.step(36, "And user enters \"MC38 guide plan\" into \"Name\" input", () => enterInto(page, "MC38 guide plan", el("\"Name\" input")), undefined, "The plan saved: the study will be followed against it");
    await session.step(38, "And user clicks on \"Save\" button in \"Save plan\" dialog", () => clickOn(page, el("\"Save\" button in \"Save plan\" dialog")), undefined, "The plan saved: the study will be followed against it");
    await session.step(40, "When user clicks on overview step and the stage opens", () => clickOpensStage(page, el("overview step")), undefined, "The study under way: its plan and its measurements so far");
    await session.step(42, "And user clicks on open a study strip", () => clickOn(page, el("open a study strip")), undefined, "The study under way: its plan and its measurements so far");
    await session.step(44, "And user clicks on saved plans link", () => clickOn(page, el("saved plans link")), undefined, "The study under way: its plan and its measurements so far");
    await session.step(46, "And user clicks on guide plan menu item", () => clickOn(page, el("guide plan menu item")), undefined, "The study under way: its plan and its measurements so far");
    await session.step(48, "And user uploads \"../files/example/mc38-ld4172-anti-pd1.csv\" through second file zone", () => uploadThrough(page, "../files/example/mc38-ld4172-anti-pd1.csv", el("second file zone")), undefined, "The study under way: its plan and its measurements so far");
    await session.step(50, "And user clicks on still running choice", () => clickOn(page, el("still running choice")), undefined, "Monitor: the running study against what its plan predicted");
    await session.step(52, "And user clicks on \"Open in Monitor\" button in \"Open a study\" dialog and the stage opens", () => clickOpensStage(page, el("\"Open in Monitor\" button in \"Open a study\" dialog")), undefined, "Monitor: the running study against what its plan predicted");
    await session.step(54, "Then \"Needs attention\" section should be visible", () => shouldBe(page, el("\"Needs attention\" section"), "visible"), undefined, "Monitor: the running study against what its plan predicted");
    await session.step(56, "When user clicks on \"Measurements\" tab", () => clickOn(page, el("\"Measurements\" tab")), undefined, "Each animal against the range expected next; any of them flagged for review, with its reason");
    await session.step(58, "And user picks \"Flag this animal…\" from the context menu of first sheet row", () => pickFromContextMenu(page, "Flag this animal…", el("first sheet row")), undefined, "Each animal against the range expected next; any of them flagged for review, with its reason");
    await session.step(60, "And user enters \"Ulcerated tumor: check at the next session\" into \"Reason\" input", () => enterInto(page, "Ulcerated tumor: check at the next session", el("\"Reason\" input")), undefined, "Each animal against the range expected next; any of them flagged for review, with its reason");
    await session.step(62, "And user clicks on \"OK\" button and PMX Studio settles", () => clickSettles(page, el("\"OK\" button")), undefined, "Each animal against the range expected next; any of them flagged for review, with its reason");
    await session.step(64, "When user clicks on \"Welfare\" tab and PMX Studio settles", () => clickSettles(page, el("\"Welfare\" tab")), undefined, "Welfare against the humane endpoints, and every decision logged with who, when and why");
    await session.step(66, "When user clicks on \"Decisions\" tab and PMX Studio settles", () => clickSettles(page, el("\"Decisions\" tab")), undefined, "Welfare against the humane endpoints, and every decision logged with who, when and why");
    await session.step(68, "When user clicks on analyze header step and the stage opens", () => clickOpensStage(page, el("analyze header step")), undefined, "Analyze: the finished study, every arm against the control");
    await session.step(70, "And user points out experimental groups section", () => pointOut(page, el("experimental groups section")), undefined, "Analyze: the finished study, every arm against the control");
    await session.step(72, "And user scrolls the page to km chart", () => scrollTo(page, el("km chart")), undefined, "Time to the endpoint, projected past the last measurement");
    await session.step(74, "And user drags the slider of km slider to 24", () => dragSliderTo(page, el("km slider"), 24), undefined, "Time to the endpoint, projected past the last measurement");
    await session.step(76, "And user scrolls the page to each animal tab", () => scrollTo(page, el("each animal tab")), undefined, "Each animal: its response and its tumor over time");
    await session.step(78, "When user clicks on each animal tab and PMX Studio settles", () => clickSettles(page, el("each animal tab")), undefined, "Each animal: its response and its tumor over time");
    await session.step(80, "And user points out waterfall chart", () => pointOut(page, el("waterfall chart")), undefined, "Each animal: its response and its tumor over time");
    await session.step(82, "And user scrolls the page to twin section", () => scrollTo(page, el("twin section")), undefined, "Each animal: its response and its tumor over time");
    await session.step(84, "And user clicks on third animal row and PMX Studio settles", () => clickSettles(page, el("third animal row")), undefined, "Its digital twin: the same tumor untreated; the gap is the effect of the treatment");
    await session.step(86, "And user points out twin effect line", () => pointOut(page, el("twin effect line")), undefined, "Its digital twin: the same tumor untreated; the gap is the effect of the treatment");
    await session.step(88, "And user scrolls the page to tolerability tab", () => scrollTo(page, el("tolerability tab")), undefined, "Tolerability, and the study against its plan");
    await session.step(90, "When user clicks on \"Tolerability\" tab and PMX Studio settles", () => clickSettles(page, el("\"Tolerability\" tab")), undefined, "Tolerability, and the study against its plan");
    await session.step(92, "When user clicks on \"Against the plan\" tab and PMX Studio settles", () => clickSettles(page, el("\"Against the plan\" tab")), undefined, "Tolerability, and the study against its plan");
    await session.step(94, "When user clicks on overview step and the stage opens", () => clickOpensStage(page, el("overview step")), undefined, "The whole portfolio, from design to analysis");
  });
});
