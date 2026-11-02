/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/browse-files.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.browse]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, collapse, downloadContains, downloadThrough, expand, fileDownloaded, followingShouldBe, isExpanded, shouldBe, shouldContainText, watchDownloads} from '@datagrok-libraries/bdd/bindings/common/steps';
import {counterMatchesFolder} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {browsePanelOpen, contextPanelOpen, contextPanelShows, refreshBrowse, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {fileOnServer} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noBalloons, noErrors, openContextMenu, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Files section of the Browse tree", () => {
  const session = feature(test, "features/browse/browse-files.feature", import.meta.url);
  test("The Files section lists its file shares", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(35, "Given user is logged in", () => loggedIn(page));
    await session.step(36, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(37, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(40, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["Files---App-Data tree node inside browse tree"],["Files---Demo tree node inside browse tree"]]), [["Files---App-Data tree node inside browse tree"],["Files---Demo tree node inside browse tree"]]);
    await session.step(43, "And no errors should have been logged", () => noErrors(page));
    await session.step(44, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A folder opens as a folder view of its own", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(35, "Given user is logged in", () => loggedIn(page));
    await session.step(36, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(37, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(47, "When user clicks on Files---Demo tree node inside browse tree", () => clickOn(page, el("Files---Demo tree node inside browse tree")));
    await session.step(48, "Then the \"Demo\" view should be current", () => viewIsCurrent(page, "Demo"));
    await session.step(50, "And gallery should contain text \"chem\"", () => shouldContainText(page, el("gallery"), "chem"));
    await session.step(51, "And the gallery counter should show as many items as the \"System:DemoFiles/\" folder holds on the server", () => counterMatchesFolder(page, "System:DemoFiles/"));
    await session.step(52, "And no errors should have been logged", () => noErrors(page));
    await session.step(53, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A tabular file opens as a preview with its rows", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(35, "Given user is logged in", () => loggedIn(page));
    await session.step(36, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(37, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(56, "Given Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(57, "When user clicks on Files---Demo---demog.csv tree node inside browse tree", () => clickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
    await session.step(58, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
    await session.step(59, "And grid should show 5850 rows", () => showsRows(page, el("grid"), 5850));
    await session.step(60, "And no errors should have been logged", () => noErrors(page));
    await session.step(61, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Download from a file's menu hands over the file", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(35, "Given user is logged in", () => loggedIn(page));
    await session.step(36, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(37, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(64, "Given user watches downloads", () => watchDownloads(page));
    await session.step(65, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(66, "When user opens the context menu of Files---Demo---demog.csv tree node inside browse tree", () => openContextMenu(page, el("Files---Demo---demog.csv tree node inside browse tree")));
    await session.step(67, "And user downloads a file through \"Download\" menu item in context menu", () => downloadThrough(page, el("\"Download\" menu item in context menu")));
    await session.step(68, "Then a file \"demog.csv\" should have been downloaded", () => fileDownloaded(page, "demog.csv"));
    await session.step(69, "And the downloaded file \"demog.csv\" should contain text \"USUBJID\"", () => downloadContains(page, "demog.csv", "USUBJID"));
    await session.step(70, "And no errors should have been logged", () => noErrors(page));
    await session.step(71, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Files written on the server show in the tree after Refresh, and files of no readable format are clicked with no error", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(35, "Given user is logged in", () => loggedIn(page));
    await session.step(36, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(37, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(74, "Given a file \"System:AppData/UsageAnalysis/bdd-browse-refresh.txt\" with text \"written by a feature\" is on the server", () => fileOnServer(page, "System:AppData/UsageAnalysis/bdd-browse-refresh.txt", "written by a feature"));
    await session.step(75, "And a file \"System:AppData/UsageAnalysis/bdd-edge-odd-{time}.qqq\" with text \"an unreadable format\" is on the server", () => fileOnServer(page, session.text("System:AppData/UsageAnalysis/bdd-edge-odd-{time}.qqq"), "an unreadable format"));
    await session.step(76, "And a file \"System:AppData/UsageAnalysis/bdd-edge-empty-{time}.csv\" with text \"\" is on the server", () => fileOnServer(page, session.text("System:AppData/UsageAnalysis/bdd-edge-empty-{time}.csv"), ""));
    await session.step(77, "And a file \"System:AppData/UsageAnalysis/bdd-edge-empty-{time}.txt\" with text \"\" is on the server", () => fileOnServer(page, session.text("System:AppData/UsageAnalysis/bdd-edge-empty-{time}.txt"), ""));
    await session.step(78, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(79, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(80, "When user refreshes the browse tree", () => refreshBrowse(page));
    await session.step(83, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(84, "And user expands Files---App-Data---UsageAnalysis tree node inside browse tree", () => expand(page, el("Files---App-Data---UsageAnalysis tree node inside browse tree")));
    await session.step(85, "Then Files---App-Data---UsageAnalysis---bdd-browse-refresh.txt tree node inside browse tree should be visible", () => shouldBe(page, el("Files---App-Data---UsageAnalysis---bdd-browse-refresh.txt tree node inside browse tree"), "visible"));
    await session.step(86, "And Files---App-Data---UsageAnalysis---bdd-edge-odd-{time}.qqq tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Files---App-Data---UsageAnalysis---bdd-edge-odd-{time}.qqq tree node inside browse tree")), "visible"));
    await session.step(87, "Given user watches downloads", () => watchDownloads(page));
    await session.step(88, "When user opens the context menu of Files---App-Data---UsageAnalysis---bdd-edge-odd-{time}.qqq tree node inside browse tree", () => openContextMenu(page, el(session.text("Files---App-Data---UsageAnalysis---bdd-edge-odd-{time}.qqq tree node inside browse tree"))));
    await session.step(89, "And user downloads a file through \"Download\" menu item in context menu", () => downloadThrough(page, el("\"Download\" menu item in context menu")));
    await session.step(90, "Then a file \"bdd-edge-odd-{time}.qqq\" should have been downloaded", () => fileDownloaded(page, session.text("bdd-edge-odd-{time}.qqq")));
    await session.step(91, "And the downloaded file \"bdd-edge-odd-{time}.qqq\" should contain text \"an unreadable format\"", () => downloadContains(page, session.text("bdd-edge-odd-{time}.qqq"), "an unreadable format"));
    await session.step(92, "When user clicks on Files---App-Data---UsageAnalysis---bdd-edge-odd-{time}.qqq tree node inside browse tree", () => clickOn(page, el(session.text("Files---App-Data---UsageAnalysis---bdd-edge-odd-{time}.qqq tree node inside browse tree"))));
    await session.step(93, "Then the context panel should show \"bdd-edge-odd-{time}.qqq\"", () => contextPanelShows(page, session.text("bdd-edge-odd-{time}.qqq")));
    await session.step(94, "When user clicks on Files---App-Data---UsageAnalysis---bdd-edge-empty-{time}.csv tree node inside browse tree", () => clickOn(page, el(session.text("Files---App-Data---UsageAnalysis---bdd-edge-empty-{time}.csv tree node inside browse tree"))));
    await session.step(95, "Then the context panel should show \"bdd-edge-empty-{time}.csv\"", () => contextPanelShows(page, session.text("bdd-edge-empty-{time}.csv")));
    await session.step(96, "When user clicks on Files---App-Data---UsageAnalysis---bdd-edge-empty-{time}.txt tree node inside browse tree", () => clickOn(page, el(session.text("Files---App-Data---UsageAnalysis---bdd-edge-empty-{time}.txt tree node inside browse tree"))));
    await session.step(97, "Then the context panel should show \"bdd-edge-empty-{time}.txt\"", () => contextPanelShows(page, session.text("bdd-edge-empty-{time}.txt")));
    await session.step(98, "When user collapses Files---App-Data tree node inside browse tree", () => collapse(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(99, "Then no errors should have been logged", () => noErrors(page));
    await session.step(100, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
