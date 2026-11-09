/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/landing/new-campaign.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {campaignHas, campaignOnServerClaim, templatesOnServer} from '../../bindings/landing.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnsExactly, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openApp, packageInstalled, urlShouldContain, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Hit Design new campaign: the form a template builds, its required fields, and Start", () => {
  const session = feature(test, "features/landing/new-campaign.feature", import.meta.url);
  test("Hit Design new campaign: the form a template builds, its required fields, and Start", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(19, "And the \"Hit Design\" templates are on the server:", () => templatesOnServer(page, "Hit Design", [["name","key","fields","stages","layout"],[session.text("BDD New Campaign {time}"),session.text("BDDNEW{time}"),"Owner: String required; Budget: Number","Idea, Made",""],[session.text("BDD New Layout {time}"),session.text("BDDNLY{time}"),"Layout note: String","Idea","bar chart of Stage"]]), [["name","key","fields","stages","layout"],[session.text("BDD New Campaign {time}"),session.text("BDDNEW{time}"),"Owner: String required; Budget: Number","Idea, Made",""],[session.text("BDD New Layout {time}"),session.text("BDDNLY{time}"),"Layout note: String","Idea","bar chart of Stage"]]);
    await session.step(23, "And user opens the \"Hit Design\" app", () => openApp(page, "Hit Design"));
    await run.scenario("A template's form has Campaign Name, the template's fields and START", async () => {
      await session.step(26, "When user selects \"BDD New Campaign {time}\" in Template input", () => selectIn(page, session.text("BDD New Campaign {time}"), el("Template input")));
      await session.step(27, "Then \"Campaign Name\" input should be visible", () => shouldBe(page, el("\"Campaign Name\" input"), "visible"));
      await session.step(28, "And \"Owner\" input should be visible", () => shouldBe(page, el("\"Owner\" input"), "visible"));
      await session.step(29, "And \"Budget\" input should be visible", () => shouldBe(page, el("\"Budget\" input"), "visible"));
      await session.step(30, "And Start button should be visible", () => shouldBe(page, el("Start button"), "visible"));
    });
    await run.scenario("While a required field is empty, START is disabled", async () => {
      await session.step(33, "When user closes all views", () => closeAllViews(page));
      await session.step(34, "And user opens the \"Hit Design\" app", () => openApp(page, "Hit Design"));
      await session.step(35, "And user selects \"BDD New Campaign {time}\" in Template input", () => selectIn(page, session.text("BDD New Campaign {time}"), el("Template input")));
      await session.step(36, "Then \"Owner\" input should be visible", () => shouldBe(page, el("\"Owner\" input"), "visible"));
      await session.step(37, "And Start button should be disabled", () => shouldBe(page, el("Start button"), "disabled"));
      await session.step(38, "When user types \"Ann\" into \"Owner\" input", () => typeInto(page, "Ann", el("\"Owner\" input")));
      await session.step(39, "Then Start button should be enabled", () => shouldBe(page, el("Start button"), "enabled"));
    });
    await run.scenario("START makes the campaign under its code and name and opens it on one row with the sketcher", async () => {
      await session.step(42, "When user closes all views", () => closeAllViews(page));
      await session.step(43, "And user opens the \"Hit Design\" app", () => openApp(page, "Hit Design"));
      await session.step(44, "And user selects \"BDD New Campaign {time}\" in Template input", () => selectIn(page, session.text("BDD New Campaign {time}"), el("Template input")));
      await session.step(45, "And user types \"Ann\" into \"Owner\" input", () => typeInto(page, "Ann", el("\"Owner\" input")));
      await session.step(46, "And user types \"BDD New Run {time}\" into \"Campaign Name\" input", () => typeInto(page, session.text("BDD New Run {time}"), el("\"Campaign Name\" input")));
      await session.step(47, "And user clicks on Start button", () => clickOn(page, el("Start button")));
      await session.step(48, "Then the \"BDD New Run {time}\" view should be current", () => viewIsCurrent(page, session.text("BDD New Run {time}")));
      await session.step(49, "And the page address should contain \"campaignId=BDDNEW{time}-1\"", () => urlShouldContain(page, session.text("campaignId=BDDNEW{time}-1")));
      await session.step(50, "And sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
      await session.step(51, "And the table should have the columns \"Molecule, Stage, V-iD\"", () => columnsExactly(page, "Molecule, Stage, V-iD"));
      await session.step(52, "And the table should have 1 row", () => rowCount(page, 1));
      await session.step(53, "And the value of \"Stage\" column in row 1 should be \"Idea\"", () => valueInRow(page, "Stage", 1, "Idea"));
      await session.step(54, "And the \"Hit Design\" campaign \"BDDNEW{time}-1\" should be on the server", () => campaignOnServerClaim(page, "Hit Design", session.text("BDDNEW{time}-1")));
      await session.step(55, "And the \"Hit Design\" campaign \"BDDNEW{time}-1\" on the server should have the name \"BDD New Run {time}\"", () => campaignHas(page, "Hit Design", session.text("BDDNEW{time}-1"), "name", session.text("BDD New Run {time}")));
      await session.step(56, "And the \"Hit Design\" campaign \"BDDNEW{time}-1\" on the server should have the status \"In Progress\"", () => campaignHas(page, "Hit Design", session.text("BDDNEW{time}-1"), "status", "In Progress"));
      await session.step(57, "When user clicks on CANCEL button in sketcher dialog", () => clickOn(page, el("CANCEL button in sketcher dialog")));
      await session.step(58, "Then sketcher dialog should be absent", () => shouldBe(page, el("sketcher dialog"), "absent"));
    });
    await run.scenario("A new campaign from a template with a layout opens with that layout", async () => {
      await session.step(61, "When user closes all views", () => closeAllViews(page));
      await session.step(62, "And user opens the \"Hit Design\" app", () => openApp(page, "Hit Design"));
      await session.step(63, "And user selects \"BDD New Layout {time}\" in Template input", () => selectIn(page, session.text("BDD New Layout {time}"), el("Template input")));
      await session.step(64, "Then \"Layout note\" input should be visible", () => shouldBe(page, el("\"Layout note\" input"), "visible"));
      await session.step(65, "When user clicks on Start button", () => clickOn(page, el("Start button")));
      await session.step(66, "Then the \"BDDNLY{time}-1\" view should be current", () => viewIsCurrent(page, session.text("BDDNLY{time}-1")));
      await session.step(67, "And sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
      await session.step(68, "When user clicks on CANCEL button in sketcher dialog", () => clickOn(page, el("CANCEL button in sketcher dialog")));
      await session.step(69, "Then bar chart viewer should be visible", () => shouldBe(page, el("bar chart viewer"), "visible"));
    });
    run.finish();
  });
});
