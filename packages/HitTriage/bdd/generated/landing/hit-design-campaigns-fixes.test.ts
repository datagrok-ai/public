/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/landing/hit-design-campaigns-fixes.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {campaignsFromElsewhere, campaignsOnServer, cellContains, doubleClickHeader, groupCount, listed, listedInOrder, settingsDefaults, sortedBy, tableColumns, tableNotSorted, tableSorted} from '../../bindings/landing.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openApp} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Hit Design landing: the campaigns table lists every campaign (full list, order, groups, dates, refresh)", () => {
  const session = feature(test, "features/landing/hit-design-campaigns-fixes.feature", import.meta.url);
  test("Hit Design landing: the campaigns table lists every campaign (full list, order, groups, dates, refresh)", {tag: ["@journey", "@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the campaigns tables start from their defaults, put back at feature end", () => settingsDefaults(page));
    await session.step(21, "And the \"Hit Design\" campaigns are on the server:", () => campaignsOnServer(page, "Hit Design", [["code","name","status","created","modified","template","fields"],[session.text("BDDFIX{time}-1"),session.text("BDD Fixes 01 {time}"),session.text("BDD Many {time}"),"2099/3/10","20",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-2"),session.text("BDD Fixes 02 {time}"),session.text("BDD Many {time}"),"2099/3/11","19",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-3"),session.text("BDD Fixes 03 {time}"),session.text("BDD Many {time}"),"2099/3/12","18",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-4"),session.text("BDD Fixes 04 {time}"),session.text("BDD Many {time}"),"2099/3/13","17",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-5"),session.text("BDD Fixes 05 {time}"),session.text("BDD Many {time}"),"2099/3/14","16",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-6"),session.text("BDD Fixes 06 {time}"),session.text("BDD Many {time}"),"2099/3/15","15",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-7"),session.text("BDD Fixes 07 {time}"),session.text("BDD Many {time}"),"2099/3/16","14",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-8"),session.text("BDD Fixes 08 {time}"),session.text("BDD Many {time}"),"2099/3/17","13",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-9"),session.text("BDD Fixes 09 {time}"),session.text("BDD Many {time}"),"2099/3/18","12",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-10"),session.text("BDD Fixes 10 {time}"),session.text("BDD Many {time}"),"2099/3/19","11",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-11"),session.text("BDD Fixes 11 {time}"),session.text("BDD Many {time}"),"2099/3/20","10",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-12"),session.text("BDD Fixes 12 {time}"),session.text("BDD Many {time}"),"2099/3/21","9",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-13"),session.text("BDD Fixes 13 {time}"),session.text("BDD Many {time}"),"2099/3/22","8",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-14"),session.text("BDD Fixes 14 {time}"),session.text("BDD Many {time}"),"2099/3/23","7",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-15"),session.text("BDD Fixes 15 {time}"),session.text("BDD Many {time}"),"2099/3/24","6",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-16"),session.text("BDD Fixes 16 {time}"),session.text("BDD Many {time}"),"2099/3/25","5",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-17"),session.text("BDD Fixes 17 {time}"),session.text("BDD Many {time}"),"2099/3/26","4",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-18"),session.text("BDD Fixes 18 {time}"),session.text("BDD Many {time}"),"2099/3/27","3",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-19"),session.text("BDD Fixes 19 {time}"),session.text("BDD Many {time}"),"2099/3/28","2",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-20"),session.text("BDD Fixes 20 {time}"),session.text("BDD Many {time}"),"2099/3/29","1",session.text("BDD Fixes Template {time}"),"Due = date 2026-11-15"],[session.text("BDDFIX{time}-21"),session.text("BDD Fixes Rare {time}"),session.text("BDD Rare {time}"),"2099/1/10","30",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-22"),session.text("BDD Fixes Last {time}"),session.text("BDD Last {time}"),"2099/1/9","31",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-24"),session.text("BDD Fixes Legacy {time}"),session.text("BDD Legacy {time}"),"2099/1/8","32",session.text("gone: BDD Fixes Template {time}"),""]]), [["code","name","status","created","modified","template","fields"],[session.text("BDDFIX{time}-1"),session.text("BDD Fixes 01 {time}"),session.text("BDD Many {time}"),"2099/3/10","20",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-2"),session.text("BDD Fixes 02 {time}"),session.text("BDD Many {time}"),"2099/3/11","19",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-3"),session.text("BDD Fixes 03 {time}"),session.text("BDD Many {time}"),"2099/3/12","18",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-4"),session.text("BDD Fixes 04 {time}"),session.text("BDD Many {time}"),"2099/3/13","17",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-5"),session.text("BDD Fixes 05 {time}"),session.text("BDD Many {time}"),"2099/3/14","16",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-6"),session.text("BDD Fixes 06 {time}"),session.text("BDD Many {time}"),"2099/3/15","15",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-7"),session.text("BDD Fixes 07 {time}"),session.text("BDD Many {time}"),"2099/3/16","14",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-8"),session.text("BDD Fixes 08 {time}"),session.text("BDD Many {time}"),"2099/3/17","13",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-9"),session.text("BDD Fixes 09 {time}"),session.text("BDD Many {time}"),"2099/3/18","12",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-10"),session.text("BDD Fixes 10 {time}"),session.text("BDD Many {time}"),"2099/3/19","11",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-11"),session.text("BDD Fixes 11 {time}"),session.text("BDD Many {time}"),"2099/3/20","10",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-12"),session.text("BDD Fixes 12 {time}"),session.text("BDD Many {time}"),"2099/3/21","9",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-13"),session.text("BDD Fixes 13 {time}"),session.text("BDD Many {time}"),"2099/3/22","8",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-14"),session.text("BDD Fixes 14 {time}"),session.text("BDD Many {time}"),"2099/3/23","7",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-15"),session.text("BDD Fixes 15 {time}"),session.text("BDD Many {time}"),"2099/3/24","6",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-16"),session.text("BDD Fixes 16 {time}"),session.text("BDD Many {time}"),"2099/3/25","5",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-17"),session.text("BDD Fixes 17 {time}"),session.text("BDD Many {time}"),"2099/3/26","4",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-18"),session.text("BDD Fixes 18 {time}"),session.text("BDD Many {time}"),"2099/3/27","3",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-19"),session.text("BDD Fixes 19 {time}"),session.text("BDD Many {time}"),"2099/3/28","2",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-20"),session.text("BDD Fixes 20 {time}"),session.text("BDD Many {time}"),"2099/3/29","1",session.text("BDD Fixes Template {time}"),"Due = date 2026-11-15"],[session.text("BDDFIX{time}-21"),session.text("BDD Fixes Rare {time}"),session.text("BDD Rare {time}"),"2099/1/10","30",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-22"),session.text("BDD Fixes Last {time}"),session.text("BDD Last {time}"),"2099/1/9","31",session.text("BDD Fixes Template {time}"),""],[session.text("BDDFIX{time}-24"),session.text("BDD Fixes Legacy {time}"),session.text("BDD Legacy {time}"),"2099/1/8","32",session.text("gone: BDD Fixes Template {time}"),""]]);
    await session.step(46, "And the campaigns table of \"Hit Design\" is sorted by \"Created\" descending", () => sortedBy(page, "Hit Design", "Created", "descending"));
    await session.step(47, "And user opens the \"Hit Design\" app", () => openApp(page, "Hit Design"));
    await run.scenario("Every campaign is listed, and stays listed after a sort change and a refresh", async () => {
      await session.step(50, "Then the campaigns table should list \"BDD Fixes 20 {time}, BDD Fixes 01 {time}\"", () => listed(page, session.text("BDD Fixes 20 {time}, BDD Fixes 01 {time}")));
      await session.step(51, "And the campaigns table should list \"BDD Fixes Rare {time}\"", () => listed(page, session.text("BDD Fixes Rare {time}")));
      await session.step(52, "When user double-clicks the \"Created\" header of the campaigns table", () => doubleClickHeader(page, "Created"));
      await session.step(53, "Then the campaigns table should not be sorted", () => tableNotSorted(page));
      await session.step(54, "And the campaigns table should list \"BDD Fixes Rare {time}, BDD Fixes 20 {time}\"", () => listed(page, session.text("BDD Fixes Rare {time}, BDD Fixes 20 {time}")));
      await session.step(55, "When user clicks on campaigns refresh icon", () => clickOn(page, el("campaigns refresh icon")));
      await session.step(56, "Then the campaigns table should list \"BDD Fixes Rare {time}, BDD Fixes 20 {time}\"", () => listed(page, session.text("BDD Fixes Rare {time}, BDD Fixes 20 {time}")));
      await session.step(57, "When user double-clicks the \"Created\" header of the campaigns table", () => doubleClickHeader(page, "Created"));
      await session.step(58, "And user double-clicks the \"Created\" header of the campaigns table", () => doubleClickHeader(page, "Created"));
      await session.step(59, "Then the campaigns table should be sorted by \"Created\" descending", () => tableSorted(page, "Created", "descending"));
    });
    await run.scenario("Grouping covers every campaign, not only the rows shown first", async () => {
      await session.step(62, "When user clicks on campaign grouping icon", () => clickOn(page, el("campaign grouping icon")));
      await session.step(63, "And user picks \"Status\" from the open menu", () => pickFromOpenMenu(page, "Status"));
      await session.step(64, "Then \"BDD Many {time}\" campaign group should be present", () => shouldBe(page, el(session.text("\"BDD Many {time}\" campaign group")), "present"));
      await session.step(65, "And the campaigns table should list \"BDD Fixes 20 {time}\"", () => listed(page, session.text("BDD Fixes 20 {time}")));
      await session.step(66, "And \"BDD Last {time}\" campaign group should be present", () => shouldBe(page, el(session.text("\"BDD Last {time}\" campaign group")), "present"));
      await session.step(67, "When user clicks on campaign grouping icon", () => clickOn(page, el("campaign grouping icon")));
      await session.step(68, "And user picks \"None\" from the open menu", () => pickFromOpenMenu(page, "None"));
      await session.step(69, "Then the campaigns table should list \"BDD Fixes 20 {time}\"", () => listed(page, session.text("BDD Fixes 20 {time}")));
    });
    await run.scenario("Groups start expanded and show how many campaigns they hold", async () => {
      await session.step(72, "When user clicks on campaign grouping icon", () => clickOn(page, el("campaign grouping icon")));
      await session.step(73, "And user picks \"Status\" from the open menu", () => pickFromOpenMenu(page, "Status"));
      await session.step(74, "Then \"BDD Many {time}\" campaign group should be present", () => shouldBe(page, el(session.text("\"BDD Many {time}\" campaign group")), "present"));
      await session.step(75, "And \"BDD Many {time}\" campaign group should be expanded", () => shouldBe(page, el(session.text("\"BDD Many {time}\" campaign group")), "expanded"));
      await session.step(76, "And \"BDD Fixes 20 {time}\" campaign row should be visible", () => shouldBe(page, el(session.text("\"BDD Fixes 20 {time}\" campaign row")), "visible"));
      await session.step(77, "And the \"BDD Many {time}\" campaign group should count 20 campaigns", () => groupCount(page, session.text("BDD Many {time}"), 20));
      await session.step(78, "When user clicks on campaign grouping icon", () => clickOn(page, el("campaign grouping icon")));
      await session.step(79, "And user picks \"None\" from the open menu", () => pickFromOpenMenu(page, "None"));
      await session.step(80, "Then the campaigns table should list \"BDD Fixes 20 {time}\"", () => listed(page, session.text("BDD Fixes 20 {time}")));
    });
    await run.scenario("A Template group is labelled by the template's name, and holds its campaigns saved without a snapshot too", async () => {
      await session.step(83, "When user clicks on campaign grouping icon", () => clickOn(page, el("campaign grouping icon")));
      await session.step(84, "And user picks \"Template\" from the open menu", () => pickFromOpenMenu(page, "Template"));
      await session.step(85, "Then \"BDD Fixes Template {time}\" campaign group should be present", () => shouldBe(page, el(session.text("\"BDD Fixes Template {time}\" campaign group")), "present"));
      await session.step(86, "And the \"BDD Fixes Template {time}\" campaign group should count 23 campaigns", () => groupCount(page, session.text("BDD Fixes Template {time}"), 23));
      await session.step(87, "When user clicks on campaign grouping icon", () => clickOn(page, el("campaign grouping icon")));
      await session.step(88, "And user picks \"None\" from the open menu", () => pickFromOpenMenu(page, "None"));
      await session.step(89, "Then the campaigns table should list \"BDD Fixes 20 {time}\"", () => listed(page, session.text("BDD Fixes 20 {time}")));
    });
    await run.scenario("A Date campaign field shows its date", async () => {
      await session.step(92, "When user clicks on campaign grouping icon", () => clickOn(page, el("campaign grouping icon")));
      await session.step(93, "And user picks \"None\" from the open menu", () => pickFromOpenMenu(page, "None"));
      await session.step(94, "And user clicks on campaign columns icon", () => clickOn(page, el("campaign columns icon")));
      await session.step(95, "And user picks \"Campaign Fields > Due\" from the open menu", () => pickFromOpenMenu(page, "Campaign Fields > Due"));
      await session.step(96, "Then the campaigns table should have the columns \"Name, Code, Created, Author, Last Modified by, Molecules, Status, Due\"", () => tableColumns(page, "Name, Code, Created, Author, Last Modified by, Molecules, Status, Due"));
      await session.step(97, "And the \"Due\" column of the \"BDD Fixes 20 {time}\" campaign row should contain \"2026\"", () => cellContains(page, "Due", session.text("BDD Fixes 20 {time}"), "2026"));
      await session.step(98, "When user picks \"Campaign Fields > Due\" from the open menu", () => pickFromOpenMenu(page, "Campaign Fields > Due"));
      await session.step(99, "And user closes the context menu", () => closeContextMenu(page));
      await session.step(100, "Then the campaigns table should have the columns \"Name, Code, Created, Author, Last Modified by, Molecules, Status\"", () => tableColumns(page, "Name, Code, Created, Author, Last Modified by, Molecules, Status"));
    });
    await run.scenario("Refresh shows a campaign saved elsewhere after the landing opened", async () => {
      await session.step(103, "When user clicks on campaign grouping icon", () => clickOn(page, el("campaign grouping icon")));
      await session.step(104, "And user picks \"None\" from the open menu", () => pickFromOpenMenu(page, "None"));
      await session.step(105, "And the \"Hit Design\" campaigns are saved on the server from elsewhere:", () => campaignsFromElsewhere(page, "Hit Design", [["code","name","status","created","modified","template"],[session.text("BDDFIX{time}-23"),session.text("BDD Fixes Late {time}"),session.text("BDD Many {time}"),"2099/4/10","0",session.text("BDD Fixes Template {time}")]]), [["code","name","status","created","modified","template"],[session.text("BDDFIX{time}-23"),session.text("BDD Fixes Late {time}"),session.text("BDD Many {time}"),"2099/4/10","0",session.text("BDD Fixes Template {time}")]]);
      await session.step(108, "And user clicks on campaigns refresh icon", () => clickOn(page, el("campaigns refresh icon")));
      await session.step(109, "Then the campaigns table should list \"BDD Fixes Late {time}\"", () => listed(page, session.text("BDD Fixes Late {time}")));
    });
    await run.scenario("With no sort chosen, the most recently modified campaigns come first", async () => {
      await session.step(112, "When user clicks on campaign grouping icon", () => clickOn(page, el("campaign grouping icon")));
      await session.step(113, "And user picks \"None\" from the open menu", () => pickFromOpenMenu(page, "None"));
      await session.step(114, "And user double-clicks the \"Created\" header of the campaigns table", () => doubleClickHeader(page, "Created"));
      await session.step(115, "Then the campaigns table should not be sorted", () => tableNotSorted(page));
      await session.step(116, "And the campaigns table should list \"BDD Fixes 20 {time}, BDD Fixes 19 {time}, BDD Fixes 18 {time}, BDD Fixes 17 {time}, BDD Fixes 16 {time}\" in this order", () => listedInOrder(page, session.text("BDD Fixes 20 {time}, BDD Fixes 19 {time}, BDD Fixes 18 {time}, BDD Fixes 17 {time}, BDD Fixes 16 {time}")));
    });
    run.finish();
  });
});
