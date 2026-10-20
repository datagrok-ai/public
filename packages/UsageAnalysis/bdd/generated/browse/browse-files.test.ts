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
import {browsePanelOpen, refreshBrowse, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {fileOnServer} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noBalloons, noErrors, openContextMenu, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Files section of the Browse tree", () => {
  const session = feature(test, "features/browse/browse-files.feature", import.meta.url);
  test("The Files section lists its file shares", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(30, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(33, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["Files---App-Data tree node inside browse tree"],["Files---Demo tree node inside browse tree"]]), [["Files---App-Data tree node inside browse tree"],["Files---Demo tree node inside browse tree"]]);
    await session.step(36, "And no errors should have been logged", () => noErrors(page));
    await session.step(37, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A folder opens as a folder view of its own", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(30, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(40, "When user clicks on Files---Demo tree node inside browse tree", () => clickOn(page, el("Files---Demo tree node inside browse tree")));
    await session.step(41, "Then the \"Demo\" view should be current", () => viewIsCurrent(page, "Demo"));
    await session.step(43, "And gallery should contain text \"chem\"", () => shouldContainText(page, el("gallery"), "chem"));
    await session.step(44, "And the gallery counter should show as many items as the \"System:DemoFiles/\" folder holds on the server", () => counterMatchesFolder(page, "System:DemoFiles/"));
    await session.step(45, "And no errors should have been logged", () => noErrors(page));
    await session.step(46, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A tabular file opens as a preview with its rows", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(30, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(49, "Given Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(50, "When user clicks on Files---Demo---demog.csv tree node inside browse tree", () => clickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
    await session.step(51, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
    await session.step(52, "And grid should show 5850 rows", () => showsRows(page, el("grid"), 5850));
    await session.step(53, "And no errors should have been logged", () => noErrors(page));
    await session.step(54, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Download from a file's menu hands over the file", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(30, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(57, "Given user watches downloads", () => watchDownloads(page));
    await session.step(58, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(59, "When user opens the context menu of Files---Demo---demog.csv tree node inside browse tree", () => openContextMenu(page, el("Files---Demo---demog.csv tree node inside browse tree")));
    await session.step(60, "And user downloads a file through \"Download\" menu item in context menu", () => downloadThrough(page, el("\"Download\" menu item in context menu")));
    await session.step(61, "Then a file \"demog.csv\" should have been downloaded", () => fileDownloaded(page, "demog.csv"));
    await session.step(62, "And the downloaded file \"demog.csv\" should contain text \"USUBJID\"", () => downloadContains(page, "demog.csv", "USUBJID"));
    await session.step(63, "And no errors should have been logged", () => noErrors(page));
    await session.step(64, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A file written on the server shows in the tree after Refresh", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(30, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(67, "Given a file \"System:AppData/UsageAnalysis/bdd-browse-refresh.txt\" with text \"written by a feature\" is on the server", () => fileOnServer(page, "System:AppData/UsageAnalysis/bdd-browse-refresh.txt", "written by a feature"));
    await session.step(68, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(69, "When user refreshes the browse tree", () => refreshBrowse(page));
    await session.step(72, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(73, "And user expands Files---App-Data---UsageAnalysis tree node inside browse tree", () => expand(page, el("Files---App-Data---UsageAnalysis tree node inside browse tree")));
    await session.step(74, "Then Files---App-Data---UsageAnalysis---bdd-browse-refresh.txt tree node inside browse tree should be visible", () => shouldBe(page, el("Files---App-Data---UsageAnalysis---bdd-browse-refresh.txt tree node inside browse tree"), "visible"));
    await session.step(75, "When user collapses Files---App-Data tree node inside browse tree", () => collapse(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(76, "Then no errors should have been logged", () => noErrors(page));
    await session.step(77, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
