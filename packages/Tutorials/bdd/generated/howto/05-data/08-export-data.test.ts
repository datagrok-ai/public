/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/05-data/08-export-data.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
import '../../../bindings/howto-databases.js';
import '../../../bindings/howto-filters.js';
import '../../../bindings/howto-ml.js';
import '../../../bindings/howto-molecules.js';
import '../../../bindings/howto-more-viewers.js';
import '../../../bindings/howto-scripting.js';
import '../../../bindings/howto-sequences.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, downloadThrough, downloadedContains, downloadedNotContains, lastDownloadCount, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPasses} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Data", () => {
  const session = feature(test, "features/howto/05-data/08-export-data.feature", import.meta.url);
  test("Export data", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(17, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(18, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(20, "When user clicks on Export icon in toolbar", () => clickOn(page, el("Export icon in toolbar")), undefined, "Export: download the table as a file");
    await session.step(22, "And user downloads a file through \"As CSV\" text in toolbar", () => downloadThrough(page, el("\"As CSV\" text in toolbar")), undefined, "Pick a format; the file goes to your downloads");
    await session.step(23, "Then the downloaded file should contain 434 occurrences of \",RA,\"", () => lastDownloadCount(page, 434, ",RA,"));
    await session.step(24, "And the downloaded file should contain \"Psoriasis\"", () => downloadedContains(page, "Psoriasis"));
    await session.step(26, "When user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")), undefined, "Filter the rows you want to keep");
    await session.step(27, "And user clicks on the \"category RA of DIS_POP\" area of filter panel", () => clickArea(page, "category RA of DIS_POP", el("filter panel")));
    await session.step(28, "Then 434 rows should pass the filter", () => filterPasses(page, 434));
    await session.step(30, "When user clicks on Export icon in toolbar", () => clickOn(page, el("Export icon in toolbar")), undefined, "CSV with options: what goes into the file");
    await session.step(31, "And user clicks on \"As CSV (options)...\" text in toolbar", () => clickOn(page, el("\"As CSV (options)...\" text in toolbar")));
    await session.step(32, "Then \"Save as CSV\" dialog should be visible", () => shouldBe(page, el("\"Save as CSV\" dialog"), "visible"));
    await session.step(34, "When user checks \"Filtered rows only\" input in \"Save as CSV\" dialog", () => check(page, el("\"Filtered rows only\" input in \"Save as CSV\" dialog")), undefined, "Filtered rows only: just the 434 rows on screen");
    await session.step(35, "And user downloads a file through OK button in \"Save as CSV\" dialog", () => downloadThrough(page, el("OK button in \"Save as CSV\" dialog")));
    await session.step(36, "Then the downloaded file should contain 434 occurrences of \",RA,\"", () => lastDownloadCount(page, 434, ",RA,"));
    await session.step(37, "And the downloaded file should not contain \"Psoriasis\"", () => downloadedNotContains(page, "Psoriasis"));
    await session.step(38, "And no errors should have been logged", () => noErrors(page));
  });
});
