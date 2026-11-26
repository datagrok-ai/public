/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/10-ml-and-statistics/02-pca.feature
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
import {gridColorsText} from '../../../bindings/howto-ml.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe, shouldContainText, shouldHaveValue, typeInto, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {colorCategorical} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pointerAway, propertyShouldBe, viewerAdded} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInColumnSelector, toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · ML", () => {
  const session = feature(test, "features/howto/10-ml-and-statistics/02-pca.feature", import.meta.url);
  test("Reduce many columns with PCA", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And user opens winequality dataset", () => openDataset(page, ds("winequality")));
    await session.step(18, "And user colors \"type\" column categorically:", () => colorCategorical(page, "type", [["red","#b2182b"],["white","#d9a520"]]), [["red","#b2182b"],["white","#d9a520"]]);
    await session.step(21, "And the grid paints the colors of \"type\" column on its text", () => gridColorsText(page, "type"));
    await session.step(22, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(23, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(24, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(25, "When the video says \"Eleven lab measurements per wine: too many to plot at once\" at the bottom", () => videoSaysAt(page, "Eleven lab measurements per wine: too many to plot at once", "bottom"));
    await session.step(27, "And user picks \"ML > Analyze > PCA...\" from the top menu", () => pickFromTopMenu(page, "ML > Analyze > PCA..."), undefined, "PCA is under ML > Analyze");
    await session.step(28, "Then \"PCA\" dialog should be visible", () => shouldBe(page, el("\"PCA\" dialog"), "visible"));
    await session.step(30, "When user clicks on editor of Features input in \"PCA\" dialog", () => clickOn(page, el("editor of Features input in \"PCA\" dialog")), undefined, "Features: the numeric columns to combine");
    await session.step(31, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(32, "When user clicks on \"All\" link in \"Select columns...\" dialog", () => clickOn(page, el("\"All\" link in \"Select columns...\" dialog")));
    await session.step(34, "And user types \"quality\" into Search input in \"Select columns...\" dialog", () => typeInto(page, "quality", el("Search input in \"Select columns...\" dialog")), undefined, "Leave out quality: it is a score, not a measurement");
    await session.step(35, "And user toggles the \"quality\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "quality", el("\"Select columns...\" dialog")));
    await session.step(36, "Then \"Select columns...\" dialog should contain text \"11 checked\"", () => shouldContainText(page, el("\"Select columns...\" dialog"), "11 checked"));
    await session.step(37, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(38, "Then editor of Features input in \"PCA\" dialog should contain text \"(11)\"", () => shouldContainText(page, el("editor of Features input in \"PCA\" dialog"), "(11)"));
    await session.step(40, "When user hovers over Components input in \"PCA\" dialog", () => hoverOver(page, el("Components input in \"PCA\" dialog")), undefined, "Components: how many new columns to make");
    await session.step(41, "Then Components input in \"PCA\" dialog should have value \"2\"", () => shouldHaveValue(page, el("Components input in \"PCA\" dialog"), "2"));
    await session.step(43, "When user clicks on OK button in \"PCA\" dialog", () => clickOn(page, el("OK button in \"PCA\" dialog")), undefined, "OK adds them to the table as PC1 and PC2");
    await session.step(44, "Then \"PCA\" dialog should be hidden", () => shouldBe(page, el("\"PCA\" dialog"), "hidden"));
    await session.step(45, "And the table should have a column \"PC1\"", () => hasColumn(page, "PC1"));
    await session.step(46, "And the table should have a column \"PC2\"", () => hasColumn(page, "PC2"));
    await session.step(47, "And \"PC1\" column should have no missing values", () => columnComplete(page, "PC1"));
    await session.step(49, "When user clicks on scatter plot icon on toolbox", () => clickOn(page, el("scatter plot icon on toolbox")), undefined, "Plot the two components against each other");
    await session.step(50, "Then scatter plot viewer should be added to the open tableview", () => viewerAdded(page, "scatter plot"));
    await session.step(51, "When user picks \"PC1\" in the \"x\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "PC1", "x", el("scatter plot viewer")));
    await session.step(52, "And user picks \"PC2\" in the \"y\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "PC2", "y", el("scatter plot viewer")));
    await session.step(53, "Then \"X\" property of scatter plot viewer should be \"PC1\"", () => propertyShouldBe(page, "X", el("scatter plot viewer"), "PC1"));
    await session.step(54, "And \"Y\" property of scatter plot viewer should be \"PC2\"", () => propertyShouldBe(page, "Y", el("scatter plot viewer"), "PC2"));
    await session.step(56, "When user picks \"type\" in the \"color\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "type", "color", el("scatter plot viewer")), undefined, "Color by type: red and white wines sit apart along PC1");
    await session.step(57, "Then \"Color\" property of scatter plot viewer should be \"type\"", () => propertyShouldBe(page, "Color", el("scatter plot viewer"), "type"));
    await session.step(58, "When user moves the pointer away from scatter plot viewer", () => pointerAway(page, el("scatter plot viewer")));
    await session.step(59, "Then no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(60, "And no errors should have been logged", () => noErrors(page));
  });
});
