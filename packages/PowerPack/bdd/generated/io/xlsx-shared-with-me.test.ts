/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/io/xlsx-shared-with-me.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [powerpack.import.xlsx]
--- */
import {test} from '@playwright/test';
import '../../bindings/add-new-column.js';
import '../../bindings/enrichment.js';
import '../../bindings/formula-lines.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {expandShared, homeWidgetsAre, openShared, runningAccountKnown, sharingUserAtHand, signBackIn, signInAsSharingUser} from '../../bindings/home.js';
import {cellOfTable, fixtureInSpace, refreshBrowseTree, tableViewsOpen} from '../../bindings/io.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, isExpanded} from '@datagrok-libraries/bdd/bindings/common/steps';
import {tableColumns, tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, dialogCloses, noSpaceOnServer, pickSharingUser, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("An Excel workbook opens from Shared with me", () => {
  const session = feature(test, "features/io/xlsx-shared-with-me.feature", import.meta.url);
  test("An Excel workbook opens from Shared with me", {tag: ["@journey", "@serial", "@realizes:powerpack.import.xlsx"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the sharing user can sign in on this page", () => sharingUserAtHand(page));
    await session.step(17, "And the name of the running account is remembered", () => runningAccountKnown(page));
    await run.scenario("A space holding the workbook is shared with the sharing user", async () => {
      await session.step(20, "Given no space named \"bdd-xlsx-{time}\" is on the server", () => noSpaceOnServer(page, session.text("bdd-xlsx-{time}")));
      await session.step(21, "And the \"fixtures/xlsx-open-test.xlsx\" file of the project is in the space \"bdd-xlsx-{time}\" as \"xlsx-open-test.xlsx\"", () => fixtureInSpace(page, "fixtures/xlsx-open-test.xlsx", session.text("bdd-xlsx-{time}"), "xlsx-open-test.xlsx"));
      await session.step(22, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(23, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(24, "When user picks \"Share...\" from the context menu of \"Spaces > bdd-xlsx-{time}\" tree node inside browse tree", () => pickFromContextMenu(page, "Share...", el(session.text("\"Spaces > bdd-xlsx-{time}\" tree node inside browse tree"))));
      await session.step(25, "And user picks the sharing user in \"User, group, or email\" input in \"Share bdd-xlsx-{time}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share bdd-xlsx-{time}\" dialog"))));
      await session.step(26, "And user clicks on OK button in \"Share bdd-xlsx-{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share bdd-xlsx-{time}\" dialog"))));
      await session.step(27, "Then the \"Share bdd-xlsx-{time}\" dialog should close", () => dialogCloses(page, session.text("Share bdd-xlsx-{time}")));
      await session.step(28, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The sharing user opens the workbook from Shared with me", async () => {
      await session.step(31, "When user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(32, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(33, "And user refreshes the browse tree", () => refreshBrowseTree(page));
      await session.step(34, "And \"My stuff\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff\" tree node inside browse tree")));
      await session.step(35, "And \"My stuff > Shared with me\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff > Shared with me\" tree node inside browse tree")));
      await session.step(36, "And user expands \".\" shared by the running account", () => expandShared(page, "."));
      await session.step(37, "And user expands \"bdd-xlsx-{time}\" shared by the running account", () => expandShared(page, session.text("bdd-xlsx-{time}")));
      await session.step(38, "And user double-clicks on \"bdd-xlsx-{time} > xlsx-open-test.xlsx\" shared by the running account", () => openShared(page, session.text("bdd-xlsx-{time} > xlsx-open-test.xlsx")));
      await session.step(39, "Then the table views \"Customers, Orders, Products\" should be open", () => tableViewsOpen(page, "Customers, Orders, Products"));
      await session.step(40, "And the \"Products\" view should be current", () => viewIsCurrent(page, "Products"));
      await session.step(41, "And table \"Customers\" should have 5 rows", () => tableRows(page, "Customers", 5));
      await session.step(42, "And table \"Customers\" should have columns \"CustomerID, Name, Country, Since\"", () => tableColumns(page, "Customers", "CustomerID, Name, Country, Since"));
      await session.step(43, "And table \"Orders\" should have 6 rows", () => tableRows(page, "Orders", 6));
      await session.step(44, "And the value of \"Amount\" column in row 1 of table \"Orders\" should be \"250.50\"", () => cellOfTable(page, "Amount", 1, "Orders", "250.50"));
      await session.step(45, "And table \"Products\" should have 4 rows", () => tableRows(page, "Products", 4));
      await session.step(46, "And the value of \"InStock\" column in row 2 of table \"Products\" should be \"17\"", () => cellOfTable(page, "InStock", 2, "Products", "17"));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
      await session.step(48, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The running account signs back in", async () => {
      await session.step(51, "When user signs back in", () => signBackIn(page));
      await session.step(52, "Then the Home page should show the widgets \"Spotlight, Reports, Usage, Community\"", () => homeWidgetsAre(page, "Spotlight, Reports, Usage, Community"));
    });
    run.finish();
  });
});
