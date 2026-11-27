/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/10-ml-and-statistics/03-map-rows-with-umap.feature
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
import {clickOn, hoverOver, shouldBe, shouldContainText, shouldHaveValue, typeInto, videoSays, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {colorCategorical} from '@datagrok-libraries/bdd/bindings/platform/data';
import {taskBarFinished, watchTaskBar} from '@datagrok-libraries/bdd/bindings/platform/events';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pointerAway, propertyShouldBe} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInColumnSelector, toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · ML", () => {
  const session = feature(test, "features/howto/10-ml-and-statistics/03-map-rows-with-umap.feature", import.meta.url);
  test("Map similar rows with UMAP", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And user opens winequality dataset", () => openDataset(page, ds("winequality")));
    await session.step(19, "And user colors \"type\" column categorically:", () => colorCategorical(page, "type", [["red","#b2182b"],["white","#d9a520"]]), [["red","#b2182b"],["white","#d9a520"]]);
    await session.step(22, "And the grid paints the colors of \"type\" column on its text", () => gridColorsText(page, "type"));
    await session.step(23, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(24, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(25, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(26, "When the video says \"A map that puts wines with similar measurements close together\" at the bottom", () => videoSaysAt(page, "A map that puts wines with similar measurements close together", "bottom"));
    await session.step(28, "And user picks \"ML > Reduce Dimensionality...\" from the top menu", () => pickFromTopMenu(page, "ML > Reduce Dimensionality..."), undefined, "Reduce Dimensionality is in the ML menu");
    await session.step(29, "Then \"Reduce Dimensionality\" dialog should be visible", () => shouldBe(page, el("\"Reduce Dimensionality\" dialog"), "visible"));
    await session.step(31, "When user clicks on Columns input in \"Reduce Dimensionality\" dialog", () => clickOn(page, el("Columns input in \"Reduce Dimensionality\" dialog")), undefined, "Columns: what to compare the rows by");
    await session.step(32, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(33, "When user clicks on \"All\" link in \"Select columns...\" dialog", () => clickOn(page, el("\"All\" link in \"Select columns...\" dialog")));
    await session.step(35, "And user types \"quality\" into Search input in \"Select columns...\" dialog", () => typeInto(page, "quality", el("Search input in \"Select columns...\" dialog")), undefined, "Leave out quality and type: the map should not see them");
    await session.step(36, "And user toggles the \"quality\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "quality", el("\"Select columns...\" dialog")));
    await session.step(37, "And user types \"type\" into Search input in \"Select columns...\" dialog", () => typeInto(page, "type", el("Search input in \"Select columns...\" dialog")));
    await session.step(38, "And user toggles the \"type\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "type", el("\"Select columns...\" dialog")));
    await session.step(39, "Then \"Select columns...\" dialog should contain text \"11 checked\"", () => shouldContainText(page, el("\"Select columns...\" dialog"), "11 checked"));
    await session.step(40, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(42, "And user hovers over Method input in \"Reduce Dimensionality\" dialog", () => hoverOver(page, el("Method input in \"Reduce Dimensionality\" dialog")), undefined, "Method: UMAP, or t-SNE");
    await session.step(43, "Then Method input in \"Reduce Dimensionality\" dialog should have value \"UMAP\"", () => shouldHaveValue(page, el("Method input in \"Reduce Dimensionality\" dialog"), "UMAP"));
    await session.step(44, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(46, "When user clicks on OK button in \"Reduce Dimensionality\" dialog", () => clickOn(page, el("OK button in \"Reduce Dimensionality\" dialog")), undefined, "OK computes the map in the browser and plots it");
    await session.step(47, "Then the task bar should have finished \"dimensionality reduction\"", () => taskBarFinished(page, "dimensionality reduction"));
    await session.step(48, "And a new column \"Cluster (DBSCAN)\" should have been added", () => newColumnNamed(page, "Cluster (DBSCAN)"));
    await session.step(49, "And \"xColumnName\" property of scatter plot viewer should be \"Embed_X_1\"", () => propertyShouldBe(page, "xColumnName", el("scatter plot viewer"), "Embed_X_1"));
    await session.step(50, "And \"colorColumnName\" property of scatter plot viewer should be \"Cluster (DBSCAN)\"", () => propertyShouldBe(page, "colorColumnName", el("scatter plot viewer"), "Cluster (DBSCAN)"));
    await session.step(51, "When the video says \"One dot per wine: close dots have similar measurements\"", () => videoSays(page, "One dot per wine: close dots have similar measurements"));
    await session.step(53, "And user picks \"type\" in the \"color\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "type", "color", el("scatter plot viewer")), undefined, "Color by type: red and white wines land in different parts of the map");
    await session.step(54, "Then \"Color\" property of scatter plot viewer should be \"type\"", () => propertyShouldBe(page, "Color", el("scatter plot viewer"), "type"));
    await session.step(55, "When user moves the pointer away from scatter plot viewer", () => pointerAway(page, el("scatter plot viewer")));
    await session.step(56, "Then no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(57, "And no errors should have been logged", () => noErrors(page));
  });
});
