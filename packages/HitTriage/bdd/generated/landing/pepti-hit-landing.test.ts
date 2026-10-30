/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/landing/pepti-hit-landing.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {campaignHas, campaignOnServerClaim, campaignViewsOpen, campaignsOnServer, cellIs, settingsDefaults, sortedBy, tableColumns, templatesOnServer} from '../../bindings/landing.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldNotContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnsExactly, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {closeAllViews, currentViewShows, openApp, packageInstalled, switchView, urlShouldContain, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("PeptiHit landing and new campaign: the campaigns table, a campaign opened by its code, and START with Helm", () => {
  const session = feature(test, "features/landing/pepti-hit-landing.feature", import.meta.url);
  test("PeptiHit landing and new campaign: the campaigns table, a campaign opened by its code, and START with Helm", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the \"Bio\" package is installed", () => packageInstalled(page, "Bio"));
    await session.step(23, "And the \"Helm\" package is installed", () => packageInstalled(page, "Helm"));
    await session.step(24, "And the campaigns tables start from their defaults, put back at feature end", () => settingsDefaults(page));
    await session.step(25, "And the \"PeptiHit\" templates are on the server:", () => templatesOnServer(page, "PeptiHit", [["name","key","fields","stages","layout"],[session.text("BDD PH Start {time}"),session.text("BDDPHS{time}"),"Start note: String","Design, Test",""],[session.text("BDD PH Named {time}"),session.text("BDDPHN{time}"),"Named note: String","Design",""],[session.text("BDD PH Fresh {time}"),session.text("BDDPHF{time}"),"Fresh note: String","Design",""],[session.text("BDD PH After {time}"),session.text("BDDPHA{time}"),"After note: String","Design",""],[session.text("BDD PH Layout {time}"),session.text("BDDPHY{time}"),"Layout note: String","Design","scatter plot of MolWt"]]), [["name","key","fields","stages","layout"],[session.text("BDD PH Start {time}"),session.text("BDDPHS{time}"),"Start note: String","Design, Test",""],[session.text("BDD PH Named {time}"),session.text("BDDPHN{time}"),"Named note: String","Design",""],[session.text("BDD PH Fresh {time}"),session.text("BDDPHF{time}"),"Fresh note: String","Design",""],[session.text("BDD PH After {time}"),session.text("BDDPHA{time}"),"After note: String","Design",""],[session.text("BDD PH Layout {time}"),session.text("BDDPHY{time}"),"Layout note: String","Design","scatter plot of MolWt"]]);
    await session.step(32, "And the \"PeptiHit\" campaigns are on the server:", () => campaignsOnServer(page, "PeptiHit", [["code","name","status","created","molecules","author","can edit"],[session.text("BDDPHL{time}-1"),session.text("BDD PH Alpha {time}"),session.text("BDD Synthesis {time}"),"2099/1/14","3","me","everyone"],[session.text("BDDPHL{time}-2"),session.text("BDD PH Beta {time}"),"In Progress","2099/1/13","2","me","everyone"],[session.text("BDDPHL{time}-3"),session.text("BDD PH Locked {time}"),"In Progress","2099/1/12","2","sharing user","author"]]), [["code","name","status","created","molecules","author","can edit"],[session.text("BDDPHL{time}-1"),session.text("BDD PH Alpha {time}"),session.text("BDD Synthesis {time}"),"2099/1/14","3","me","everyone"],[session.text("BDDPHL{time}-2"),session.text("BDD PH Beta {time}"),"In Progress","2099/1/13","2","me","everyone"],[session.text("BDDPHL{time}-3"),session.text("BDD PH Locked {time}"),"In Progress","2099/1/12","2","sharing user","author"]]);
    await session.step(37, "And the campaigns table of \"PeptiHit\" is sorted by \"Created\" descending", () => sortedBy(page, "PeptiHit", "Created", "descending"));
    await session.step(38, "And user opens the \"PeptiHit\" app", () => openApp(page, "PeptiHit"));
    await run.scenario("The PeptiHit landing lists the campaigns as Hit Design does, under its own hero", async () => {
      await session.step(41, "Then the \"PeptiHit\" view should be current", () => viewIsCurrent(page, "PeptiHit"));
      await session.step(42, "And the current view should show the text \"Sketch Helm sequences in the spreadsheet\"", () => currentViewShows(page, "Sketch Helm sequences in the spreadsheet"));
      await session.step(43, "And the campaigns table should have the columns \"Name, Code, Created, Author, Last Modified by, Molecules, Status\"", () => tableColumns(page, "Name, Code, Created, Author, Last Modified by, Molecules, Status"));
      await session.step(44, "And the \"Status\" column of the \"BDD PH Alpha {time}\" campaign row should be \"BDD Synthesis {time}\"", () => cellIs(page, "Status", session.text("BDD PH Alpha {time}"), session.text("BDD Synthesis {time}")));
      await session.step(45, "And permissions icon of \"BDD PH Alpha {time}\" campaign row should be visible", () => shouldBe(page, el(session.text("permissions icon of \"BDD PH Alpha {time}\" campaign row")), "visible"));
      await session.step(46, "And permissions icon of \"BDD PH Locked {time}\" campaign row should be hidden", () => shouldBe(page, el(session.text("permissions icon of \"BDD PH Locked {time}\" campaign row")), "hidden"));
    });
    await run.scenario("A campaign's name opens it with the Helm column first, and the crumb comes back", async () => {
      await session.step(49, "When user clicks on \"BDD PH Alpha {time}\" link", () => clickOn(page, el(session.text("\"BDD PH Alpha {time}\" link"))));
      await session.step(50, "Then the page address should contain \"campaignId=BDDPHL{time}-1\"", () => urlShouldContain(page, session.text("campaignId=BDDPHL{time}-1")));
      await session.step(51, "And 1 campaign view should be open", () => campaignViewsOpen(page, 1));
      await session.step(52, "And the table should have the columns \"Helm, Molecule, Stage, V-iD\"", () => columnsExactly(page, "Helm, Molecule, Stage, V-iD"));
      await session.step(53, "And \"Helm\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "Helm", "Macromolecule"));
      await session.step(54, "When user clicks on \"PeptiHit\" link", () => clickOn(page, el("\"PeptiHit\" link")));
      await session.step(55, "Then the \"PeptiHit\" view should be current", () => viewIsCurrent(page, "PeptiHit"));
    });
    await run.scenario("START makes a campaign with Helm, Molecule, Stage and V-iD and opens the HELM editor on it", async () => {
      await session.step(58, "When user closes all views", () => closeAllViews(page));
      await session.step(59, "And user opens the \"PeptiHit\" app", () => openApp(page, "PeptiHit"));
      await session.step(60, "And user selects \"BDD PH Start {time}\" in Template input", () => selectIn(page, session.text("BDD PH Start {time}"), el("Template input")));
      await session.step(61, "Then \"Start note\" input should be visible", () => shouldBe(page, el("\"Start note\" input"), "visible"));
      await session.step(62, "When user clicks on Start button", () => clickOn(page, el("Start button")));
      await session.step(63, "Then HELM editor dialog should be visible", () => shouldBe(page, el("HELM editor dialog"), "visible"));
      await session.step(64, "And the \"BDDPHS{time}-1\" view should be current", () => viewIsCurrent(page, session.text("BDDPHS{time}-1")));
      await session.step(65, "And the page address should contain \"campaignId=BDDPHS{time}-1\"", () => urlShouldContain(page, session.text("campaignId=BDDPHS{time}-1")));
      await session.step(66, "And the table should have the columns \"Helm, Molecule, Stage, V-iD\"", () => columnsExactly(page, "Helm, Molecule, Stage, V-iD"));
      await session.step(67, "And the value of \"Stage\" column in row 1 should be \"Design\"", () => valueInRow(page, "Stage", 1, "Design"));
      await session.step(68, "And \"Helm\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "Helm", "Macromolecule"));
      await session.step(69, "And the \"PeptiHit\" campaign \"BDDPHS{time}-1\" should be on the server", () => campaignOnServerClaim(page, "PeptiHit", session.text("BDDPHS{time}-1")));
      await session.step(70, "When user clicks on CANCEL button in HELM editor dialog", () => clickOn(page, el("CANCEL button in HELM editor dialog")));
      await session.step(71, "Then HELM editor dialog should be absent", () => shouldBe(page, el("HELM editor dialog"), "absent"));
    });
    await run.scenario("The Campaign Name typed before START names the new campaign", async () => {
      await session.step(74, "When user closes all views", () => closeAllViews(page));
      await session.step(75, "And user opens the \"PeptiHit\" app", () => openApp(page, "PeptiHit"));
      await session.step(76, "And user selects \"BDD PH Named {time}\" in Template input", () => selectIn(page, session.text("BDD PH Named {time}"), el("Template input")));
      await session.step(77, "Then \"Named note\" input should be visible", () => shouldBe(page, el("\"Named note\" input"), "visible"));
      await session.step(78, "When user types \"BDD PH Named Run {time}\" into \"Campaign Name\" input", () => typeInto(page, session.text("BDD PH Named Run {time}"), el("\"Campaign Name\" input")));
      await session.step(79, "And user clicks on Start button", () => clickOn(page, el("Start button")));
      await session.step(80, "Then HELM editor dialog should be visible", () => shouldBe(page, el("HELM editor dialog"), "visible"));
      await session.step(81, "And the \"PeptiHit\" campaign \"BDDPHN{time}-1\" should be on the server", () => campaignOnServerClaim(page, "PeptiHit", session.text("BDDPHN{time}-1")));
      await session.step(82, "And the \"PeptiHit\" campaign \"BDDPHN{time}-1\" on the server should have the name \"BDD PH Named Run {time}\"", () => campaignHas(page, "PeptiHit", session.text("BDDPHN{time}-1"), "name", session.text("BDD PH Named Run {time}")));
      await session.step(83, "When user clicks on CANCEL button in HELM editor dialog", () => clickOn(page, el("CANCEL button in HELM editor dialog")));
    });
    await run.scenario("A new campaign does not take over the campaign opened before it", async () => {
      await session.step(86, "When user closes all views", () => closeAllViews(page));
      await session.step(87, "And user opens the \"PeptiHit\" app", () => openApp(page, "PeptiHit"));
      await session.step(88, "And user clicks on \"BDD PH Alpha {time}\" link", () => clickOn(page, el(session.text("\"BDD PH Alpha {time}\" link"))));
      await session.step(89, "Then the page address should contain \"campaignId=BDDPHL{time}-1\"", () => urlShouldContain(page, session.text("campaignId=BDDPHL{time}-1")));
      await session.step(90, "When user clicks on \"PeptiHit\" link", () => clickOn(page, el("\"PeptiHit\" link")));
      await session.step(91, "And user selects \"BDD PH Fresh {time}\" in Template input", () => selectIn(page, session.text("BDD PH Fresh {time}"), el("Template input")));
      await session.step(92, "Then \"Fresh note\" input should be visible", () => shouldBe(page, el("\"Fresh note\" input"), "visible"));
      await session.step(93, "When user clicks on Start button", () => clickOn(page, el("Start button")));
      await session.step(94, "Then HELM editor dialog should be visible", () => shouldBe(page, el("HELM editor dialog"), "visible"));
      await session.step(95, "And the \"PeptiHit\" campaign \"BDDPHF{time}-1\" should be on the server", () => campaignOnServerClaim(page, "PeptiHit", session.text("BDDPHF{time}-1")));
      await session.step(96, "And the \"PeptiHit\" campaign \"BDDPHF{time}-1\" on the server should have the status \"In Progress\"", () => campaignHas(page, "PeptiHit", session.text("BDDPHF{time}-1"), "status", "In Progress"));
      await session.step(97, "When user clicks on CANCEL button in HELM editor dialog", () => clickOn(page, el("CANCEL button in HELM editor dialog")));
    });
    await run.scenario("After a campaign the user may only view, START still makes and opens the new campaign", async () => {
      await session.step(100, "When user closes all views", () => closeAllViews(page));
      await session.step(101, "And user opens the \"PeptiHit\" app", () => openApp(page, "PeptiHit"));
      await session.step(102, "And user clicks on \"BDD PH Locked {time}\" link", () => clickOn(page, el(session.text("\"BDD PH Locked {time}\" link"))));
      await session.step(103, "Then the page address should contain \"campaignId=BDDPHL{time}-3\"", () => urlShouldContain(page, session.text("campaignId=BDDPHL{time}-3")));
      await session.step(104, "When user clicks on \"PeptiHit\" link", () => clickOn(page, el("\"PeptiHit\" link")));
      await session.step(105, "And user selects \"BDD PH After {time}\" in Template input", () => selectIn(page, session.text("BDD PH After {time}"), el("Template input")));
      await session.step(106, "Then \"After note\" input should be visible", () => shouldBe(page, el("\"After note\" input"), "visible"));
      await session.step(107, "When user clicks on Start button", () => clickOn(page, el("Start button")));
      await session.step(108, "Then the \"PeptiHit\" campaign \"BDDPHA{time}-1\" should be on the server", () => campaignOnServerClaim(page, "PeptiHit", session.text("BDDPHA{time}-1")));
      await session.step(109, "And the \"BDDPHA{time}-1\" view should be current", () => viewIsCurrent(page, session.text("BDDPHA{time}-1")));
      await session.step(110, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Opening a campaign from the landing closes the campaign opened before it", async () => {
      await session.step(113, "When user closes all views", () => closeAllViews(page));
      await session.step(114, "And user opens the \"PeptiHit\" app", () => openApp(page, "PeptiHit"));
      await session.step(115, "And user clicks on \"BDD PH Beta {time}\" link", () => clickOn(page, el(session.text("\"BDD PH Beta {time}\" link"))));
      await session.step(116, "Then the page address should contain \"campaignId=BDDPHL{time}-2\"", () => urlShouldContain(page, session.text("campaignId=BDDPHL{time}-2")));
      await session.step(117, "When user switches to the \"PeptiHit\" view", () => switchView(page, "PeptiHit"));
      await session.step(118, "And user clicks on \"BDD PH Locked {time}\" link", () => clickOn(page, el(session.text("\"BDD PH Locked {time}\" link"))));
      await session.step(119, "Then the page address should contain \"campaignId=BDDPHL{time}-3\"", () => urlShouldContain(page, session.text("campaignId=BDDPHL{time}-3")));
      await session.step(120, "And 1 campaign view should be open", () => campaignViewsOpen(page, 1));
    });
    await run.scenario("A new campaign from a template with a layout opens with viewers that work", async () => {
      await session.step(123, "When user closes all views", () => closeAllViews(page));
      await session.step(124, "And user opens the \"PeptiHit\" app", () => openApp(page, "PeptiHit"));
      await session.step(125, "And user selects \"BDD PH Layout {time}\" in Template input", () => selectIn(page, session.text("BDD PH Layout {time}"), el("Template input")));
      await session.step(126, "Then \"Layout note\" input should be visible", () => shouldBe(page, el("\"Layout note\" input"), "visible"));
      await session.step(127, "When user clicks on Start button", () => clickOn(page, el("Start button")));
      await session.step(128, "Then HELM editor dialog should be visible", () => shouldBe(page, el("HELM editor dialog"), "visible"));
      await session.step(129, "When user clicks on CANCEL button in HELM editor dialog", () => clickOn(page, el("CANCEL button in HELM editor dialog")));
      await session.step(130, "Then the \"BDDPHY{time}-1\" view should be current", () => viewIsCurrent(page, session.text("BDDPHY{time}-1")));
      await session.step(131, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(132, "And scatter plot viewer should not contain text \"does not exist\"", () => shouldNotContainText(page, el("scatter plot viewer"), "does not exist"));
    });
    run.finish();
  });
});
