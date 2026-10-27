/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/file-preview-edge-cases.feature
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
import {clickOn, downloadContains, downloadThrough, fileDownloaded, isExpanded, shouldBe, watchDownloads} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelOpen, contextPanelShows, refreshBrowse} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {fileOnServer} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noBalloons, noErrors, openContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("File preview edge cases: a file of a format the platform does not read, an empty file", () => {
  const session = feature(test, "features/browse/file-preview-edge-cases.feature", import.meta.url);
  test("A file of an unsupported format is clicked with no error and shows its details", {tag: ["@browse", "@edge-cases", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(16, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(17, "And a file \"System:AppData/UsageAnalysis/bdd-edge-odd-{time}.qqq\" with text \"an unreadable format\" is on the server", () => fileOnServer(page, session.text("System:AppData/UsageAnalysis/bdd-edge-odd-{time}.qqq"), "an unreadable format"));
    await session.step(18, "And a file \"System:AppData/UsageAnalysis/bdd-edge-empty-{time}.csv\" with text \"\" is on the server", () => fileOnServer(page, session.text("System:AppData/UsageAnalysis/bdd-edge-empty-{time}.csv"), ""));
    await session.step(19, "And a file \"System:AppData/UsageAnalysis/bdd-edge-empty-{time}.txt\" with text \"\" is on the server", () => fileOnServer(page, session.text("System:AppData/UsageAnalysis/bdd-edge-empty-{time}.txt"), ""));
    await session.step(20, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(21, "When user refreshes the browse tree", () => refreshBrowse(page));
    await session.step(22, "Given Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(23, "And Files---App-Data---UsageAnalysis tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---UsageAnalysis tree node inside browse tree")));
    await session.step(24, "Then Files---App-Data---UsageAnalysis---bdd-edge-odd-{time}.qqq tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Files---App-Data---UsageAnalysis---bdd-edge-odd-{time}.qqq tree node inside browse tree")), "visible"));
    await session.step(27, "When user clicks on Files---App-Data---UsageAnalysis---bdd-edge-odd-{time}.qqq tree node inside browse tree", () => clickOn(page, el(session.text("Files---App-Data---UsageAnalysis---bdd-edge-odd-{time}.qqq tree node inside browse tree"))));
    await session.step(28, "Then the context panel should show \"bdd-edge-odd-{time}.qqq\"", () => contextPanelShows(page, session.text("bdd-edge-odd-{time}.qqq")));
    await session.step(29, "And no errors should have been logged", () => noErrors(page));
    await session.step(30, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A file of an unsupported format is downloaded from its menu", {tag: ["@browse", "@edge-cases", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(16, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(17, "And a file \"System:AppData/UsageAnalysis/bdd-edge-odd-{time}.qqq\" with text \"an unreadable format\" is on the server", () => fileOnServer(page, session.text("System:AppData/UsageAnalysis/bdd-edge-odd-{time}.qqq"), "an unreadable format"));
    await session.step(18, "And a file \"System:AppData/UsageAnalysis/bdd-edge-empty-{time}.csv\" with text \"\" is on the server", () => fileOnServer(page, session.text("System:AppData/UsageAnalysis/bdd-edge-empty-{time}.csv"), ""));
    await session.step(19, "And a file \"System:AppData/UsageAnalysis/bdd-edge-empty-{time}.txt\" with text \"\" is on the server", () => fileOnServer(page, session.text("System:AppData/UsageAnalysis/bdd-edge-empty-{time}.txt"), ""));
    await session.step(20, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(21, "When user refreshes the browse tree", () => refreshBrowse(page));
    await session.step(22, "Given Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(23, "And Files---App-Data---UsageAnalysis tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---UsageAnalysis tree node inside browse tree")));
    await session.step(24, "Then Files---App-Data---UsageAnalysis---bdd-edge-odd-{time}.qqq tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Files---App-Data---UsageAnalysis---bdd-edge-odd-{time}.qqq tree node inside browse tree")), "visible"));
    await session.step(33, "Given user watches downloads", () => watchDownloads(page));
    await session.step(34, "When user opens the context menu of Files---App-Data---UsageAnalysis---bdd-edge-odd-{time}.qqq tree node inside browse tree", () => openContextMenu(page, el(session.text("Files---App-Data---UsageAnalysis---bdd-edge-odd-{time}.qqq tree node inside browse tree"))));
    await session.step(35, "And user downloads a file through \"Download\" menu item in context menu", () => downloadThrough(page, el("\"Download\" menu item in context menu")));
    await session.step(36, "Then a file \"bdd-edge-odd-{time}.qqq\" should have been downloaded", () => fileDownloaded(page, session.text("bdd-edge-odd-{time}.qqq")));
    await session.step(37, "And the downloaded file \"bdd-edge-odd-{time}.qqq\" should contain text \"an unreadable format\"", () => downloadContains(page, session.text("bdd-edge-odd-{time}.qqq"), "an unreadable format"));
    await session.step(38, "And no errors should have been logged", () => noErrors(page));
    await session.step(39, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("An empty csv file is previewed with no error [ext=csv]", {tag: ["@browse", "@edge-cases", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(16, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(17, "And a file \"System:AppData/UsageAnalysis/bdd-edge-odd-{time}.qqq\" with text \"an unreadable format\" is on the server", () => fileOnServer(page, session.text("System:AppData/UsageAnalysis/bdd-edge-odd-{time}.qqq"), "an unreadable format"));
    await session.step(18, "And a file \"System:AppData/UsageAnalysis/bdd-edge-empty-{time}.csv\" with text \"\" is on the server", () => fileOnServer(page, session.text("System:AppData/UsageAnalysis/bdd-edge-empty-{time}.csv"), ""));
    await session.step(19, "And a file \"System:AppData/UsageAnalysis/bdd-edge-empty-{time}.txt\" with text \"\" is on the server", () => fileOnServer(page, session.text("System:AppData/UsageAnalysis/bdd-edge-empty-{time}.txt"), ""));
    await session.step(20, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(21, "When user refreshes the browse tree", () => refreshBrowse(page));
    await session.step(22, "Given Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(23, "And Files---App-Data---UsageAnalysis tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---UsageAnalysis tree node inside browse tree")));
    await session.step(24, "Then Files---App-Data---UsageAnalysis---bdd-edge-odd-{time}.qqq tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Files---App-Data---UsageAnalysis---bdd-edge-odd-{time}.qqq tree node inside browse tree")), "visible"));
    await session.step(42, "When user clicks on Files---App-Data---UsageAnalysis---bdd-edge-empty-{time}.csv tree node inside browse tree", () => clickOn(page, el(session.text("Files---App-Data---UsageAnalysis---bdd-edge-empty-{time}.csv tree node inside browse tree"))));
    await session.step(43, "Then the context panel should show \"bdd-edge-empty-{time}.csv\"", () => contextPanelShows(page, session.text("bdd-edge-empty-{time}.csv")));
    await session.step(44, "And no errors should have been logged", () => noErrors(page));
    await session.step(45, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("An empty txt file is previewed with no error [ext=txt]", {tag: ["@browse", "@edge-cases", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(16, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(17, "And a file \"System:AppData/UsageAnalysis/bdd-edge-odd-{time}.qqq\" with text \"an unreadable format\" is on the server", () => fileOnServer(page, session.text("System:AppData/UsageAnalysis/bdd-edge-odd-{time}.qqq"), "an unreadable format"));
    await session.step(18, "And a file \"System:AppData/UsageAnalysis/bdd-edge-empty-{time}.csv\" with text \"\" is on the server", () => fileOnServer(page, session.text("System:AppData/UsageAnalysis/bdd-edge-empty-{time}.csv"), ""));
    await session.step(19, "And a file \"System:AppData/UsageAnalysis/bdd-edge-empty-{time}.txt\" with text \"\" is on the server", () => fileOnServer(page, session.text("System:AppData/UsageAnalysis/bdd-edge-empty-{time}.txt"), ""));
    await session.step(20, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(21, "When user refreshes the browse tree", () => refreshBrowse(page));
    await session.step(22, "Given Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(23, "And Files---App-Data---UsageAnalysis tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---UsageAnalysis tree node inside browse tree")));
    await session.step(24, "Then Files---App-Data---UsageAnalysis---bdd-edge-odd-{time}.qqq tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Files---App-Data---UsageAnalysis---bdd-edge-odd-{time}.qqq tree node inside browse tree")), "visible"));
    await session.step(42, "When user clicks on Files---App-Data---UsageAnalysis---bdd-edge-empty-{time}.txt tree node inside browse tree", () => clickOn(page, el(session.text("Files---App-Data---UsageAnalysis---bdd-edge-empty-{time}.txt tree node inside browse tree"))));
    await session.step(43, "Then the context panel should show \"bdd-edge-empty-{time}.txt\"", () => contextPanelShows(page, session.text("bdd-edge-empty-{time}.txt")));
    await session.step(44, "And no errors should have been logged", () => noErrors(page));
    await session.step(45, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
