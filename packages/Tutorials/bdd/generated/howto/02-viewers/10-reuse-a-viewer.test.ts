/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/02-viewers/10-reuse-a-viewer.feature
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
import {videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hasArea, hasNoArea, noErrors, pickFromAreaContextMenu, propertyShouldBe, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Viewers", () => {
  const session = feature(test, "features/howto/02-viewers/10-reuse-a-viewer.feature", import.meta.url);
  test("Copy a viewer or its look", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(18, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(19, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(20, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","HEIGHT"],["Y","WEIGHT"],["Marker Type","circle border"],["Show X Histogram","true"],["Show Y Histogram","true"]]), [["X","HEIGHT"],["Y","WEIGHT"],["Marker Type","circle border"],["Show X Histogram","true"],["Show Y Histogram","true"]]);
    await session.step(26, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","AGE"],["Y","HEIGHT"]]), [["X","AGE"],["Y","HEIGHT"]]);
    await session.step(29, "Then the open tableview should have 2 scatter plot viewers", () => viewerCount(page, 2, "scatter plot"));
    await session.step(30, "And second scatter plot viewer should not have an \"x histogram\" area", () => hasNoArea(page, el("second scatter plot viewer"), "x histogram"));
    await session.step(31, "When the video says \"Reuse the look of a viewer you have set up\"", () => videoSays(page, "Reuse the look of a viewer you have set up"));
    await session.step(33, "And user picks \"Pick Up / Apply > Pick Up\" from the context menu of the \"view\" area of first scatter plot viewer", () => pickFromAreaContextMenu(page, "Pick Up / Apply > Pick Up", "view", el("first scatter plot viewer")), undefined, "Right-click > Pick Up / Apply > Pick Up remembers its settings");
    await session.step(35, "And user picks \"Pick Up / Apply > Apply Style Settings\" from the context menu of the \"view\" area of second scatter plot viewer", () => pickFromAreaContextMenu(page, "Pick Up / Apply > Apply Style Settings", "view", el("second scatter plot viewer")), undefined, "Apply Style Settings copies the look, not the columns");
    await session.step(36, "Then second scatter plot viewer should have an \"x histogram\" area", () => hasArea(page, el("second scatter plot viewer"), "x histogram"));
    await session.step(37, "And second scatter plot viewer should have a \"y histogram\" area", () => hasArea(page, el("second scatter plot viewer"), "y histogram"));
    await session.step(38, "And \"Marker Type\" property of second scatter plot viewer should be \"circle border\"", () => propertyShouldBe(page, "Marker Type", el("second scatter plot viewer"), "circle border"));
    await session.step(39, "And \"X\" property of second scatter plot viewer should be \"AGE\"", () => propertyShouldBe(page, "X", el("second scatter plot viewer"), "AGE"));
    await session.step(41, "When user picks \"General > Clone\" from the context menu of the \"view\" area of second scatter plot viewer", () => pickFromAreaContextMenu(page, "General > Clone", "view", el("second scatter plot viewer")), undefined, "General > Clone adds a copy of a viewer");
    await session.step(42, "Then the open tableview should have 3 scatter plot viewers", () => viewerCount(page, 3, "scatter plot"));
    await session.step(43, "And \"X\" property of third scatter plot viewer should be \"AGE\"", () => propertyShouldBe(page, "X", el("third scatter plot viewer"), "AGE"));
    await session.step(44, "And third scatter plot viewer should have an \"x histogram\" area", () => hasArea(page, el("third scatter plot viewer"), "x histogram"));
    await session.step(45, "And no errors should have been logged", () => noErrors(page));
  });
});
