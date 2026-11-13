/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/landing/hit-triage-campaigns.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {campaignNotOnServer, campaignsOnServer, cellIs, confirmDialog, doubleClickHeader, menuItemChecked, settingsDefaults, sortedBy, tableColumns, tableNotSorted, tableSorted} from '../../bindings/landing.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, collapse, expand, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentViewShows, dialogCloses, openApp} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, noBalloons, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Hit Triage landing: the campaigns table lists, sorts, groups and acts on campaigns", () => {
  const session = feature(test, "features/landing/hit-triage-campaigns.feature", import.meta.url);
  test("Hit Triage landing: the campaigns table lists, sorts, groups and acts on campaigns", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the campaigns tables start from their defaults, put back at feature end", () => settingsDefaults(page));
    await session.step(18, "And the \"Hit Triage\" campaigns are on the server:", () => campaignsOnServer(page, "Hit Triage", [["code","name","status","created","molecules","author","can edit"],[session.text("BDDHTL{time}-1"),session.text("BDD HT Alpha {time}"),session.text("BDD Ready {time}"),"2099/1/13","4","me","everyone"],[session.text("BDDHTL{time}-2"),session.text("BDD HT Beta {time}"),session.text("BDD Hold {time}"),"2099/1/12","9","me","everyone"],[session.text("BDDHTL{time}-3"),session.text("BDD HT Gamma {time}"),session.text("BDD Ready {time}"),"2099/1/11","6","sharing user","author"]]), [["code","name","status","created","molecules","author","can edit"],[session.text("BDDHTL{time}-1"),session.text("BDD HT Alpha {time}"),session.text("BDD Ready {time}"),"2099/1/13","4","me","everyone"],[session.text("BDDHTL{time}-2"),session.text("BDD HT Beta {time}"),session.text("BDD Hold {time}"),"2099/1/12","9","me","everyone"],[session.text("BDDHTL{time}-3"),session.text("BDD HT Gamma {time}"),session.text("BDD Ready {time}"),"2099/1/11","6","sharing user","author"]]);
    await session.step(23, "And the campaigns table of \"Hit Triage\" is sorted by \"Created\" descending", () => sortedBy(page, "Hit Triage", "Created", "descending"));
    await session.step(24, "And user opens the \"Hit Triage\" app", () => openApp(page, "Hit Triage"));
    await run.scenario("The table lists the campaigns with Hit Triage's default columns", async () => {
      await session.step(27, "Then the current view should show the text \"Continue Campaign\"", () => currentViewShows(page, "Continue Campaign"));
      await session.step(28, "And the campaigns table should have the columns \"Name, Created, Author, Total, Selected, Status\"", () => tableColumns(page, "Name, Created, Author, Total, Selected, Status"));
      await session.step(29, "And the \"Total\" column of the \"BDD HT Beta {time}\" campaign row should be \"9\"", () => cellIs(page, "Total", session.text("BDD HT Beta {time}"), "9"));
      await session.step(30, "And the \"Selected\" column of the \"BDD HT Beta {time}\" campaign row should be \"9\"", () => cellIs(page, "Selected", session.text("BDD HT Beta {time}"), "9"));
      await session.step(31, "And the \"Status\" column of the \"BDD HT Gamma {time}\" campaign row should be \"BDD Ready {time}\"", () => cellIs(page, "Status", session.text("BDD HT Gamma {time}"), session.text("BDD Ready {time}")));
      await session.step(32, "And name link of \"BDD HT Alpha {time}\" campaign row should be visible", () => shouldBe(page, el(session.text("name link of \"BDD HT Alpha {time}\" campaign row")), "visible"));
    });
    await run.scenario("The eye menu adds a column at the end and takes it away", async () => {
      await session.step(35, "When user clicks on campaign columns icon", () => clickOn(page, el("campaign columns icon")));
      await session.step(36, "Then the open menu should show \"Code\" as unchecked", () => menuItemChecked(page, "Code", "unchecked"));
      await session.step(37, "And the open menu should show \"Total\" as checked", () => menuItemChecked(page, "Total", "checked"));
      await session.step(38, "When user picks \"Code\" from the open menu", () => pickFromOpenMenu(page, "Code"));
      await session.step(39, "Then the campaigns table should have the columns \"Name, Created, Author, Total, Selected, Status, Code\"", () => tableColumns(page, "Name, Created, Author, Total, Selected, Status, Code"));
      await session.step(40, "And the \"Code\" column of the \"BDD HT Alpha {time}\" campaign row should be \"BDDHTL{time}-1\"", () => cellIs(page, "Code", session.text("BDD HT Alpha {time}"), session.text("BDDHTL{time}-1")));
      await session.step(41, "When user picks \"Code\" from the open menu", () => pickFromOpenMenu(page, "Code"));
      await session.step(42, "Then the campaigns table should have the columns \"Name, Created, Author, Total, Selected, Status\"", () => tableColumns(page, "Name, Created, Author, Total, Selected, Status"));
      await session.step(43, "When user closes the context menu", () => closeContextMenu(page));
    });
    await run.scenario("A double click on a header sorts up, then down, then not at all", async () => {
      await session.step(46, "Then the campaigns table should be sorted by \"Created\" descending", () => tableSorted(page, "Created", "descending"));
      await session.step(47, "When user double-clicks the \"Total\" header of the campaigns table", () => doubleClickHeader(page, "Total"));
      await session.step(48, "Then the campaigns table should be sorted by \"Total\" ascending", () => tableSorted(page, "Total", "ascending"));
      await session.step(49, "When user double-clicks the \"Total\" header of the campaigns table", () => doubleClickHeader(page, "Total"));
      await session.step(50, "Then the campaigns table should be sorted by \"Total\" descending", () => tableSorted(page, "Total", "descending"));
      await session.step(51, "When user double-clicks the \"Total\" header of the campaigns table", () => doubleClickHeader(page, "Total"));
      await session.step(52, "Then the campaigns table should not be sorted", () => tableNotSorted(page));
      await session.step(53, "When user double-clicks the \"Created\" header of the campaigns table", () => doubleClickHeader(page, "Created"));
      await session.step(54, "And user double-clicks the \"Created\" header of the campaigns table", () => doubleClickHeader(page, "Created"));
      await session.step(55, "Then the campaigns table should be sorted by \"Created\" descending", () => tableSorted(page, "Created", "descending"));
    });
    await run.scenario("Grouping by Status puts a header over each status's rows, which expands and collapses them", async () => {
      await session.step(58, "When user clicks on campaign grouping icon", () => clickOn(page, el("campaign grouping icon")));
      await session.step(59, "And user picks \"Status\" from the open menu", () => pickFromOpenMenu(page, "Status"));
      await session.step(60, "Then \"BDD Ready {time}\" campaign group should be visible", () => shouldBe(page, el(session.text("\"BDD Ready {time}\" campaign group")), "visible"));
      await session.step(61, "And \"BDD Hold {time}\" campaign group should be visible", () => shouldBe(page, el(session.text("\"BDD Hold {time}\" campaign group")), "visible"));
      await session.step(62, "When user expands \"BDD Ready {time}\" campaign group", () => expand(page, el(session.text("\"BDD Ready {time}\" campaign group"))));
      await session.step(63, "Then \"BDD HT Alpha {time}\" campaign row should be visible", () => shouldBe(page, el(session.text("\"BDD HT Alpha {time}\" campaign row")), "visible"));
      await session.step(64, "And \"BDD HT Gamma {time}\" campaign row should be visible", () => shouldBe(page, el(session.text("\"BDD HT Gamma {time}\" campaign row")), "visible"));
      await session.step(65, "When user collapses \"BDD Ready {time}\" campaign group", () => collapse(page, el(session.text("\"BDD Ready {time}\" campaign group"))));
      await session.step(66, "Then \"BDD HT Alpha {time}\" campaign row should be hidden", () => shouldBe(page, el(session.text("\"BDD HT Alpha {time}\" campaign row")), "hidden"));
      await session.step(67, "When user clicks on campaign grouping icon", () => clickOn(page, el("campaign grouping icon")));
      await session.step(68, "And user picks \"None\" from the open menu", () => pickFromOpenMenu(page, "None"));
      await session.step(69, "Then \"BDD HT Alpha {time}\" campaign row should be visible", () => shouldBe(page, el(session.text("\"BDD HT Alpha {time}\" campaign row")), "visible"));
      await session.step(70, "And \"BDD Ready {time}\" campaign group should be absent", () => shouldBe(page, el(session.text("\"BDD Ready {time}\" campaign group")), "absent"));
    });
    await run.scenario("Only the campaigns the user may edit show the share and delete icons", async () => {
      await session.step(73, "Then permissions icon of \"BDD HT Alpha {time}\" campaign row should be visible", () => shouldBe(page, el(session.text("permissions icon of \"BDD HT Alpha {time}\" campaign row")), "visible"));
      await session.step(74, "And delete icon of \"BDD HT Alpha {time}\" campaign row should be visible", () => shouldBe(page, el(session.text("delete icon of \"BDD HT Alpha {time}\" campaign row")), "visible"));
      await session.step(75, "And permissions icon of \"BDD HT Gamma {time}\" campaign row should be hidden", () => shouldBe(page, el(session.text("permissions icon of \"BDD HT Gamma {time}\" campaign row")), "hidden"));
      await session.step(76, "And delete icon of \"BDD HT Gamma {time}\" campaign row should be hidden", () => shouldBe(page, el(session.text("delete icon of \"BDD HT Gamma {time}\" campaign row")), "hidden"));
    });
    await run.scenario("The share icon opens Edit Permissions, and CANCEL closes it", async () => {
      await session.step(79, "When user clicks on permissions icon of \"BDD HT Alpha {time}\" campaign row", () => clickOn(page, el(session.text("permissions icon of \"BDD HT Alpha {time}\" campaign row"))));
      await session.step(80, "Then \"Edit Permissions\" dialog should be visible", () => shouldBe(page, el("\"Edit Permissions\" dialog"), "visible"));
      await session.step(81, "And \"Edit Permissions\" dialog should contain text \"Can View\"", () => shouldContainText(page, el("\"Edit Permissions\" dialog"), "Can View"));
      await session.step(82, "And \"Edit Permissions\" dialog should contain text \"Can Edit\"", () => shouldContainText(page, el("\"Edit Permissions\" dialog"), "Can Edit"));
      await session.step(83, "When user clicks on CANCEL button in \"Edit Permissions\" dialog", () => clickOn(page, el("CANCEL button in \"Edit Permissions\" dialog")));
      await session.step(84, "Then the \"Edit Permissions\" dialog should close", () => dialogCloses(page, "Edit Permissions"));
    });
    await run.scenario("Delete asks for a confirmation naming the campaign, and OK removes it", async () => {
      await session.step(87, "When user clicks on delete icon of \"BDD HT Beta {time}\" campaign row", () => clickOn(page, el(session.text("delete icon of \"BDD HT Beta {time}\" campaign row"))));
      await session.step(88, "Then \"Delete campaign\" dialog should contain text \"Are you sure you want to delete campaign\"", () => shouldContainText(page, el("\"Delete campaign\" dialog"), "Are you sure you want to delete campaign"));
      await session.step(89, "And \"Delete campaign\" dialog should contain text \"BDD HT Beta {time}\"", () => shouldContainText(page, el("\"Delete campaign\" dialog"), session.text("BDD HT Beta {time}")));
      await session.step(90, "When user confirms the \"Delete campaign\" dialog", () => confirmDialog(page, el("the \"Delete campaign\" dialog")));
      await session.step(91, "Then the \"Hit Triage\" campaign \"BDDHTL{time}-2\" should not be on the server", () => campaignNotOnServer(page, "Hit Triage", session.text("BDDHTL{time}-2")));
      await session.step(92, "And \"BDD HT Alpha {time}\" campaign row should be visible", () => shouldBe(page, el(session.text("\"BDD HT Alpha {time}\" campaign row")), "visible"));
      await session.step(93, "And \"BDD HT Beta {time}\" campaign row should be absent", () => shouldBe(page, el(session.text("\"BDD HT Beta {time}\" campaign row")), "absent"));
    });
    await run.scenario("The refresh icon rebuilds the table in the order it had", async () => {
      await session.step(96, "When user clicks on campaigns refresh icon", () => clickOn(page, el("campaigns refresh icon")));
      await session.step(97, "Then the campaigns table should be sorted by \"Created\" descending", () => tableSorted(page, "Created", "descending"));
      await session.step(98, "And \"BDD HT Alpha {time}\" campaign row should be visible", () => shouldBe(page, el(session.text("\"BDD HT Alpha {time}\" campaign row")), "visible"));
      await session.step(99, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
