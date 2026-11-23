/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/render/reaction-export.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.reactions.export]
--- */
import {test} from '@playwright/test';
import '../../bindings/crux.js';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '../../bindings/reaction-editor.js';
import '../../bindings/reactions-app-text.js';
import '../../bindings/sketcher-app.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {engineGate} from '../../bindings/engine.js';
import {choose, coloursPerPreset, exportFromRibbon, openAppOnDefault, panelOffers, previewAtPrintSize, previewSmaller, previewWidth, save, savedPng} from '../../bindings/reaction-export.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, expand, fileDownloaded, watchDownloads} from '@datagrok-libraries/bdd/bindings/common/steps';
import {contextPanelOpen, contextPanelShowsCurrentCell, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors, pickFromAreaContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Export panel previews a reaction at print size and saves what it previews", () => {
  const session = feature(test, "features/render/reaction-export.feature", import.meta.url);
  test("The Export panel previews a reaction at print size and saves what it previews", {tag: ["@journey", "@realizes:chem.reactions.export", "@eval-b7"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await run.scenario("From the Reactions app's ribbon, on the default route", async () => {
      await session.step(17, "Given user opens Chem's Reactions app on its default example", () => openAppOnDefault(page));
      await session.step(18, "And user watches downloads", () => watchDownloads(page));
      await session.step(19, "When user clicks Export... in the Reactions app's ribbon", () => exportFromRibbon(page));
      await session.step(20, "Then the Export panel should offer its presets, widths and formats, and no PDF", () => panelOffers(page));
      await session.step(21, "And the Export panel's centre colours should be off for \"ACS Document 1996\" and on for \"Screen\"", () => coloursPerPreset(page, "ACS Document 1996", "Screen"));
      await session.step(22, "When user chooses \"ACS Document 1996\", \"single column (8.3 cm)\" and \"PNG 300 dpi\" in the Export panel", () => choose(page, "ACS Document 1996", "single column (8.3 cm)", "PNG 300 dpi"));
      await session.step(23, "Then the Export panel's preview should be 8.3 cm wide in the \"publication\" look", () => previewWidth(page, 8.3, "publication"));
      await session.step(24, "And the Export panel's note should say how large the preview is against print size", () => previewSmaller(page));
      await session.step(25, "When user clicks Save in the Export panel", () => save(page));
      await session.step(26, "Then a file \"reaction.png\" should have been downloaded", () => fileDownloaded(page, "reaction.png"));
      await session.step(27, "And the saved PNG should be 980 px wide at 300 dpi, with no colour in it", () => savedPng(page, 980, 300));
      await session.step(28, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("From a reaction cell's menu", async () => {
      await session.step(31, "Given user opens test-reactions dataset", () => openDataset(page, ds("test-reactions")));
      await session.step(32, "And user watches downloads", () => watchDownloads(page));
      await session.step(33, "When user picks \"Export reaction...\" from the context menu of the \"cell 5 of reaction\" area of grid", () => pickFromAreaContextMenu(page, "Export reaction...", "cell 5 of reaction", el("grid")));
      await session.step(34, "And user chooses \"ACS Document 1996\", \"single column (8.3 cm)\" and \"PNG 300 dpi\" in the Export panel", () => choose(page, "ACS Document 1996", "single column (8.3 cm)", "PNG 300 dpi"));
      await session.step(35, "Then the Export panel's preview should be 8.3 cm wide in the \"publication\" look", () => previewWidth(page, 8.3, "publication"));
      await session.step(36, "And the Export panel's preview should be at print size", () => previewAtPrintSize(page));
      await session.step(37, "When user clicks Save in the Export panel", () => save(page));
      await session.step(38, "Then the saved PNG should be 980 px wide at 300 dpi, with no colour in it", () => savedPng(page, 980, 300));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("From the Reaction pane", async () => {
      await session.step(42, "Given user opens test-reactions dataset", () => openDataset(page, ds("test-reactions")));
      await session.step(43, "And the context panel is open", () => contextPanelOpen(page));
      await session.step(44, "When user clicks on the \"cell 5 of reaction\" area of grid", () => clickArea(page, "cell 5 of reaction", el("grid")));
      await session.step(45, "Then the context panel should show the current cell", () => contextPanelShowsCurrentCell(page));
      await session.step(46, "When user expands \"Reaction\" section", () => expand(page, el("\"Reaction\" section")));
      await session.step(47, "And user clicks on \"Export...\" icon inside \"Reaction\" section", () => clickOn(page, el("\"Export...\" icon inside \"Reaction\" section")));
      await session.step(48, "And user chooses \"ACS Document 1996\", \"double column (17.8 cm)\" and \"PNG 600 dpi\" in the Export panel", () => choose(page, "ACS Document 1996", "double column (17.8 cm)", "PNG 600 dpi"));
      await session.step(49, "Then the Export panel's preview should be 17.8 cm wide in the \"publication\" look", () => previewWidth(page, 17.8, "publication"));
      await session.step(50, "And the Export panel's preview should be at print size", () => previewAtPrintSize(page));
      await session.step(51, "When user clicks Save in the Export panel", () => save(page));
      await session.step(52, "Then the saved PNG should be 4205 px wide at 600 dpi, with no colour in it", () => savedPng(page, 4205, 600));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
