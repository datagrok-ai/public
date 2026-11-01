/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/landing/hit-design-campaigns.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {campaignNotOnServer, campaignsOnServer, cellIs, confirmDialog, doubleClickHeader, listedInOrder, menuItemChecked, settingsDefaults, sortedBy, tableColumns, tableNotSorted, tableSorted} from '../../bindings/landing.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, collapse, expand, shouldBe, shouldContainText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, currentViewShows, dialogCloses, openApp, urlShouldContain, urlShouldNotContain, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, menuDoesNotList, menuLists, noBalloons, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Hit Design landing: the campaigns table lists, opens, sorts, groups and acts on campaigns", () => {
  const session = feature(test, "features/landing/hit-design-campaigns.feature", import.meta.url);
  test("Hit Design landing: the campaigns table lists, opens, sorts, groups and acts on campaigns", {tag: ["@journey", "@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 12, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the campaigns tables start from their defaults, put back at feature end", () => settingsDefaults(page));
    await session.step(19, "And the \"Hit Design\" campaigns are on the server:", () => campaignsOnServer(page, "Hit Design", [["code","name","status","created","molecules","author","can edit","fields"],[session.text("BDDLND{time}-1"),session.text("BDD Landing Alpha {time}"),session.text("BDD Ready {time}"),"2099/1/14","3","me","everyone",session.text("Owner = Ann {time}")],[session.text("BDDLND{time}-2"),session.text("BDD Landing Beta {time}"),session.text("BDD Hold {time}"),"2099/1/13","7","me","everyone",session.text("Owner = Bob {time}")],[session.text("BDDLND{time}-3"),session.text("BDD Landing Gamma {time}"),session.text("BDD Ready {time}"),"2099/1/12","5","sharing user","author",session.text("Owner = Ann {time}")],[session.text("BDDLND{time}-4"),session.text("BDD Landing Delta {time}"),session.text("BDD Hold {time}"),"2099/1/11","2","me","everyone",""]]), [["code","name","status","created","molecules","author","can edit","fields"],[session.text("BDDLND{time}-1"),session.text("BDD Landing Alpha {time}"),session.text("BDD Ready {time}"),"2099/1/14","3","me","everyone",session.text("Owner = Ann {time}")],[session.text("BDDLND{time}-2"),session.text("BDD Landing Beta {time}"),session.text("BDD Hold {time}"),"2099/1/13","7","me","everyone",session.text("Owner = Bob {time}")],[session.text("BDDLND{time}-3"),session.text("BDD Landing Gamma {time}"),session.text("BDD Ready {time}"),"2099/1/12","5","sharing user","author",session.text("Owner = Ann {time}")],[session.text("BDDLND{time}-4"),session.text("BDD Landing Delta {time}"),session.text("BDD Hold {time}"),"2099/1/11","2","me","everyone",""]]);
    await session.step(25, "And the campaigns table of \"Hit Design\" is sorted by \"Created\" descending", () => sortedBy(page, "Hit Design", "Created", "descending"));
    await session.step(26, "And user opens the \"Hit Design\" app", () => openApp(page, "Hit Design"));
    await run.scenario("The table lists the campaigns with the default columns, each name a link", async () => {
      await session.step(29, "Then the current view should show the text \"Continue Campaign\"", () => currentViewShows(page, "Continue Campaign"));
      await session.step(30, "And the campaigns table should have the columns \"Name, Code, Created, Author, Last Modified by, Molecules, Status\"", () => tableColumns(page, "Name, Code, Created, Author, Last Modified by, Molecules, Status"));
      await session.step(31, "And the \"Code\" column of the \"BDD Landing Alpha {time}\" campaign row should be \"BDDLND{time}-1\"", () => cellIs(page, "Code", session.text("BDD Landing Alpha {time}"), session.text("BDDLND{time}-1")));
      await session.step(32, "And the \"Created\" column of the \"BDD Landing Alpha {time}\" campaign row should be \"2099/1/14\"", () => cellIs(page, "Created", session.text("BDD Landing Alpha {time}"), "2099/1/14"));
      await session.step(33, "And the \"Molecules\" column of the \"BDD Landing Beta {time}\" campaign row should be \"7\"", () => cellIs(page, "Molecules", session.text("BDD Landing Beta {time}"), "7"));
      await session.step(34, "And the \"Status\" column of the \"BDD Landing Gamma {time}\" campaign row should be \"BDD Ready {time}\"", () => cellIs(page, "Status", session.text("BDD Landing Gamma {time}"), session.text("BDD Ready {time}")));
      await session.step(35, "And name link of \"BDD Landing Alpha {time}\" campaign row should be visible", () => shouldBe(page, el(session.text("name link of \"BDD Landing Alpha {time}\" campaign row")), "visible"));
    });
    await run.scenario("A campaign's name opens it with its code in the address, and the crumb comes back", async () => {
      await session.step(38, "When user clicks on \"BDD Landing Alpha {time}\" link", () => clickOn(page, el(session.text("\"BDD Landing Alpha {time}\" link"))));
      await session.step(39, "Then the \"BDD Landing Alpha {time}\" view should be current", () => viewIsCurrent(page, session.text("BDD Landing Alpha {time}")));
      await session.step(40, "And the page address should contain \"campaignId=BDDLND{time}-1\"", () => urlShouldContain(page, session.text("campaignId=BDDLND{time}-1")));
      await session.step(41, "And the table should have 3 rows", () => rowCount(page, 3));
      await session.step(42, "When user clicks on \"Hit Design\" link", () => clickOn(page, el("\"Hit Design\" link")));
      await session.step(43, "Then the \"Hit Design\" view should be current", () => viewIsCurrent(page, "Hit Design"));
      await session.step(44, "And the page address should not contain \"campaignId\"", () => urlShouldNotContain(page, "campaignId"));
      await session.step(45, "And \"BDD Landing Alpha {time}\" campaign row should be visible", () => shouldBe(page, el(session.text("\"BDD Landing Alpha {time}\" campaign row")), "visible"));
    });
    await run.scenario("The eye menu shows and hides columns, a column comes back at the end, and the choice outlives a reopen", async () => {
      await session.step(48, "When user clicks on campaign columns icon", () => clickOn(page, el("campaign columns icon")));
      await session.step(49, "Then the open menu should show \"Code\" as checked", () => menuItemChecked(page, "Code", "checked"));
      await session.step(50, "And the open menu should show \"Status\" as checked", () => menuItemChecked(page, "Status", "checked"));
      await session.step(51, "And the open menu should list \"Campaign Fields > Owner\"", () => menuLists(page, "Campaign Fields > Owner"));
      await session.step(52, "When user picks \"Code\" from the open menu", () => pickFromOpenMenu(page, "Code"));
      await session.step(53, "Then the campaigns table should have the columns \"Name, Created, Author, Last Modified by, Molecules, Status\"", () => tableColumns(page, "Name, Created, Author, Last Modified by, Molecules, Status"));
      await session.step(54, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(55, "And user closes all views", () => closeAllViews(page));
      await session.step(56, "And user opens the \"Hit Design\" app", () => openApp(page, "Hit Design"));
      await session.step(57, "Then the campaigns table should have the columns \"Name, Created, Author, Last Modified by, Molecules, Status\"", () => tableColumns(page, "Name, Created, Author, Last Modified by, Molecules, Status"));
      await session.step(58, "When user clicks on campaign columns icon", () => clickOn(page, el("campaign columns icon")));
      await session.step(59, "Then the open menu should show \"Code\" as unchecked", () => menuItemChecked(page, "Code", "unchecked"));
      await session.step(60, "When user picks \"Code\" from the open menu", () => pickFromOpenMenu(page, "Code"));
      await session.step(61, "Then the campaigns table should have the columns \"Name, Created, Author, Last Modified by, Molecules, Status, Code\"", () => tableColumns(page, "Name, Created, Author, Last Modified by, Molecules, Status, Code"));
      await session.step(62, "When user picks \"Campaign Fields > Owner\" from the open menu", () => pickFromOpenMenu(page, "Campaign Fields > Owner"));
      await session.step(63, "Then the campaigns table should have the columns \"Name, Created, Author, Last Modified by, Molecules, Status, Code, Owner\"", () => tableColumns(page, "Name, Created, Author, Last Modified by, Molecules, Status, Code, Owner"));
      await session.step(64, "And the \"Owner\" column of the \"BDD Landing Alpha {time}\" campaign row should be \"Ann {time}\"", () => cellIs(page, "Owner", session.text("BDD Landing Alpha {time}"), session.text("Ann {time}")));
      await session.step(65, "And the \"Owner\" column of the \"BDD Landing Delta {time}\" campaign row should be \"\"", () => cellIs(page, "Owner", session.text("BDD Landing Delta {time}"), ""));
      await session.step(66, "When user picks \"Campaign Fields > Owner\" from the open menu", () => pickFromOpenMenu(page, "Campaign Fields > Owner"));
      await session.step(67, "Then the campaigns table should have the columns \"Name, Created, Author, Last Modified by, Molecules, Status, Code\"", () => tableColumns(page, "Name, Created, Author, Last Modified by, Molecules, Status, Code"));
      await session.step(68, "When user closes the context menu", () => closeContextMenu(page));
    });
    await run.scenario("A double click on a header sorts up, then down, then not at all", async () => {
      await session.step(71, "Then the campaigns table should be sorted by \"Created\" descending", () => tableSorted(page, "Created", "descending"));
      await session.step(72, "And the campaigns table should list \"BDD Landing Alpha {time}, BDD Landing Beta {time}, BDD Landing Gamma {time}, BDD Landing Delta {time}\" in this order", () => listedInOrder(page, session.text("BDD Landing Alpha {time}, BDD Landing Beta {time}, BDD Landing Gamma {time}, BDD Landing Delta {time}")));
      await session.step(73, "When user double-clicks the \"Molecules\" header of the campaigns table", () => doubleClickHeader(page, "Molecules"));
      await session.step(74, "Then the campaigns table should be sorted by \"Molecules\" ascending", () => tableSorted(page, "Molecules", "ascending"));
      await session.step(75, "When user double-clicks the \"Molecules\" header of the campaigns table", () => doubleClickHeader(page, "Molecules"));
      await session.step(76, "Then the campaigns table should be sorted by \"Molecules\" descending", () => tableSorted(page, "Molecules", "descending"));
      await session.step(77, "When user double-clicks the \"Molecules\" header of the campaigns table", () => doubleClickHeader(page, "Molecules"));
      await session.step(78, "Then the campaigns table should not be sorted", () => tableNotSorted(page));
      await session.step(79, "When user double-clicks the \"Created\" header of the campaigns table", () => doubleClickHeader(page, "Created"));
      await session.step(80, "Then the campaigns table should be sorted by \"Created\" ascending", () => tableSorted(page, "Created", "ascending"));
      await session.step(81, "When user double-clicks the \"Created\" header of the campaigns table", () => doubleClickHeader(page, "Created"));
      await session.step(82, "Then the campaigns table should be sorted by \"Created\" descending", () => tableSorted(page, "Created", "descending"));
    });
    await run.scenario("Grouping by Status puts a header over each status's rows, which expands and collapses them", async () => {
      await session.step(85, "When user clicks on campaign grouping icon", () => clickOn(page, el("campaign grouping icon")));
      await session.step(86, "And user picks \"Status\" from the open menu", () => pickFromOpenMenu(page, "Status"));
      await session.step(87, "Then \"BDD Ready {time}\" campaign group should be visible", () => shouldBe(page, el(session.text("\"BDD Ready {time}\" campaign group")), "visible"));
      await session.step(88, "And \"BDD Hold {time}\" campaign group should be visible", () => shouldBe(page, el(session.text("\"BDD Hold {time}\" campaign group")), "visible"));
      await session.step(89, "When user expands \"BDD Ready {time}\" campaign group", () => expand(page, el(session.text("\"BDD Ready {time}\" campaign group"))));
      await session.step(90, "Then \"BDD Landing Alpha {time}\" campaign row should be visible", () => shouldBe(page, el(session.text("\"BDD Landing Alpha {time}\" campaign row")), "visible"));
      await session.step(91, "And \"BDD Landing Gamma {time}\" campaign row should be visible", () => shouldBe(page, el(session.text("\"BDD Landing Gamma {time}\" campaign row")), "visible"));
      await session.step(92, "When user collapses \"BDD Ready {time}\" campaign group", () => collapse(page, el(session.text("\"BDD Ready {time}\" campaign group"))));
      await session.step(93, "Then \"BDD Landing Alpha {time}\" campaign row should be hidden", () => shouldBe(page, el(session.text("\"BDD Landing Alpha {time}\" campaign row")), "hidden"));
      await session.step(94, "And \"BDD Landing Gamma {time}\" campaign row should be hidden", () => shouldBe(page, el(session.text("\"BDD Landing Gamma {time}\" campaign row")), "hidden"));
    });
    await run.scenario("Grouping by a campaign field, with \"No Information\" for campaigns without it, then back to none", async () => {
      await session.step(97, "When user clicks on campaign grouping icon", () => clickOn(page, el("campaign grouping icon")));
      await session.step(98, "And user picks \"Campaign Fields > Owner\" from the open menu", () => pickFromOpenMenu(page, "Campaign Fields > Owner"));
      await session.step(99, "Then \"Ann {time}\" campaign group should be visible", () => shouldBe(page, el(session.text("\"Ann {time}\" campaign group")), "visible"));
      await session.step(100, "And \"Bob {time}\" campaign group should be visible", () => shouldBe(page, el(session.text("\"Bob {time}\" campaign group")), "visible"));
      await session.step(101, "And \"No Information\" campaign group should be visible", () => shouldBe(page, el("\"No Information\" campaign group"), "visible"));
      await session.step(102, "When user clicks on campaign grouping icon", () => clickOn(page, el("campaign grouping icon")));
      await session.step(103, "And user picks \"None\" from the open menu", () => pickFromOpenMenu(page, "None"));
      await session.step(104, "Then \"BDD Landing Alpha {time}\" campaign row should be visible", () => shouldBe(page, el(session.text("\"BDD Landing Alpha {time}\" campaign row")), "visible"));
      await session.step(105, "And \"Ann {time}\" campaign group should be absent", () => shouldBe(page, el(session.text("\"Ann {time}\" campaign group")), "absent"));
    });
    await run.scenario("The grouping menu offers a grouping under the caption of its column: Last Modified by", async () => {
      await session.step(108, "When user clicks on campaign grouping icon", () => clickOn(page, el("campaign grouping icon")));
      await session.step(109, "Then the open menu should list \"Last Modified by\"", () => menuLists(page, "Last Modified by"));
      await session.step(110, "And the open menu should not list \"Last Modified User\"", () => menuDoesNotList(page, "Last Modified User"));
      await session.step(111, "When user picks \"Last Modified by\" from the open menu", () => pickFromOpenMenu(page, "Last Modified by"));
      await session.step(112, "And user clicks on campaign grouping icon", () => clickOn(page, el("campaign grouping icon")));
      await session.step(113, "Then the open menu should show \"Last Modified by\" as checked", () => menuItemChecked(page, "Last Modified by", "checked"));
      await session.step(114, "When user picks \"None\" from the open menu", () => pickFromOpenMenu(page, "None"));
      await session.step(115, "Then \"BDD Landing Alpha {time}\" campaign row should be visible", () => shouldBe(page, el(session.text("\"BDD Landing Alpha {time}\" campaign row")), "visible"));
    });
    await run.scenario("Only the campaigns the user may edit show the permissions and delete icons", async () => {
      await session.step(118, "Then permissions icon of \"BDD Landing Alpha {time}\" campaign row should be visible", () => shouldBe(page, el(session.text("permissions icon of \"BDD Landing Alpha {time}\" campaign row")), "visible"));
      await session.step(119, "And delete icon of \"BDD Landing Alpha {time}\" campaign row should be visible", () => shouldBe(page, el(session.text("delete icon of \"BDD Landing Alpha {time}\" campaign row")), "visible"));
      await session.step(120, "And permissions icon of \"BDD Landing Gamma {time}\" campaign row should be hidden", () => shouldBe(page, el(session.text("permissions icon of \"BDD Landing Gamma {time}\" campaign row")), "hidden"));
      await session.step(121, "And delete icon of \"BDD Landing Gamma {time}\" campaign row should be hidden", () => shouldBe(page, el(session.text("delete icon of \"BDD Landing Gamma {time}\" campaign row")), "hidden"));
    });
    await run.scenario("The permissions icon opens Edit Permissions with Can View and Can Edit, and CANCEL closes it", async () => {
      await session.step(124, "When user clicks on permissions icon of \"BDD Landing Alpha {time}\" campaign row", () => clickOn(page, el(session.text("permissions icon of \"BDD Landing Alpha {time}\" campaign row"))));
      await session.step(125, "Then \"Edit Permissions\" dialog should be visible", () => shouldBe(page, el("\"Edit Permissions\" dialog"), "visible"));
      await session.step(126, "And \"Edit Permissions\" dialog should contain text \"Can View\"", () => shouldContainText(page, el("\"Edit Permissions\" dialog"), "Can View"));
      await session.step(127, "And \"Edit Permissions\" dialog should contain text \"Can Edit\"", () => shouldContainText(page, el("\"Edit Permissions\" dialog"), "Can Edit"));
      await session.step(128, "When user clicks on CANCEL button in \"Edit Permissions\" dialog", () => clickOn(page, el("CANCEL button in \"Edit Permissions\" dialog")));
      await session.step(129, "Then the \"Edit Permissions\" dialog should close", () => dialogCloses(page, "Edit Permissions"));
    });
    await run.scenario("Edit Permissions says all users can view only while no one is listed", async () => {
      await session.step(132, "When user clicks on permissions icon of \"BDD Landing Alpha {time}\" campaign row", () => clickOn(page, el(session.text("permissions icon of \"BDD Landing Alpha {time}\" campaign row"))));
      await session.step(133, "Then \"Edit Permissions\" dialog should contain text \"Can Edit includes Can View\"", () => shouldContainText(page, el("\"Edit Permissions\" dialog"), "Can Edit includes Can View"));
      await session.step(134, "And \"Edit Permissions\" dialog should not contain text \"With no one listed\"", () => shouldNotContainText(page, el("\"Edit Permissions\" dialog"), "With no one listed"));
      await session.step(135, "When user clicks on \"Remove All users\" button in \"Edit Permissions\" dialog", () => clickOn(page, el("\"Remove All users\" button in \"Edit Permissions\" dialog")));
      await session.step(136, "Then \"Edit Permissions\" dialog should contain text \"With no one listed, all users can view\"", () => shouldContainText(page, el("\"Edit Permissions\" dialog"), "With no one listed, all users can view"));
      await session.step(137, "When user clicks on CANCEL button in \"Edit Permissions\" dialog", () => clickOn(page, el("CANCEL button in \"Edit Permissions\" dialog")));
      await session.step(138, "Then the \"Edit Permissions\" dialog should close", () => dialogCloses(page, "Edit Permissions"));
    });
    await run.scenario("Delete asks for a confirmation naming the code, and confirming it removes the campaign", async () => {
      await session.step(141, "When user clicks on delete icon of \"BDD Landing Delta {time}\" campaign row", () => clickOn(page, el(session.text("delete icon of \"BDD Landing Delta {time}\" campaign row"))));
      await session.step(142, "Then \"Delete campaign\" dialog should contain text \"Are you sure you want to delete campaign\"", () => shouldContainText(page, el("\"Delete campaign\" dialog"), "Are you sure you want to delete campaign"));
      await session.step(143, "And \"Delete campaign\" dialog should contain text \"BDDLND{time}-4\"", () => shouldContainText(page, el("\"Delete campaign\" dialog"), session.text("BDDLND{time}-4")));
      await session.step(144, "When user confirms the \"Delete campaign\" dialog", () => confirmDialog(page, el("the \"Delete campaign\" dialog")));
      await session.step(145, "Then the \"Hit Design\" campaign \"BDDLND{time}-4\" should not be on the server", () => campaignNotOnServer(page, "Hit Design", session.text("BDDLND{time}-4")));
      await session.step(146, "And \"BDD Landing Alpha {time}\" campaign row should be visible", () => shouldBe(page, el(session.text("\"BDD Landing Alpha {time}\" campaign row")), "visible"));
      await session.step(147, "And \"BDD Landing Delta {time}\" campaign row should be absent", () => shouldBe(page, el(session.text("\"BDD Landing Delta {time}\" campaign row")), "absent"));
    });
    await run.scenario("The refresh icon rebuilds the table in the order it had", async () => {
      await session.step(150, "When user clicks on campaigns refresh icon", () => clickOn(page, el("campaigns refresh icon")));
      await session.step(151, "Then the campaigns table should be sorted by \"Created\" descending", () => tableSorted(page, "Created", "descending"));
      await session.step(152, "And \"BDD Landing Alpha {time}\" campaign row should be visible", () => shouldBe(page, el(session.text("\"BDD Landing Alpha {time}\" campaign row")), "visible"));
      await session.step(153, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
