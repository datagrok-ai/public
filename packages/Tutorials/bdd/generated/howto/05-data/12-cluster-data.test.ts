/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/05-data/12-cluster-data.feature
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
import {check, clickOn, enterInto, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {distinctValues, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, propertyShouldBe, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInColumnSelector, toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Data", () => {
  const session = feature(test, "features/howto/05-data/12-cluster-data.feature", import.meta.url);
  test("Cluster data", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And user opens iris dataset", () => openDataset(page, ds("iris")));
    await session.step(19, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(20, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(21, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(23, "When user picks \"ML > Cluster > Cluster Data...\" from the top menu", () => pickFromTopMenu(page, "ML > Cluster > Cluster Data..."), undefined, "Cluster Data: k-means, right in the browser");
    await session.step(24, "Then \"Cluster Data\" dialog should be visible", () => shouldBe(page, el("\"Cluster Data\" dialog"), "visible"));
    await session.step(26, "When user clicks on Features input in \"Cluster Data\" dialog", () => clickOn(page, el("Features input in \"Cluster Data\" dialog")), undefined, "Features: the columns that make rows alike");
    await session.step(27, "And user toggles the \"col 1\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "col 1", el("\"Select columns...\" dialog")));
    await session.step(28, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(30, "And user enters \"2\" into Clusters input in \"Cluster Data\" dialog", () => enterInto(page, "2", el("Clusters input in \"Cluster Data\" dialog")), undefined, "Clusters: how many groups to find");
    await session.step(32, "And user checks \"Show Scatter Plot\" input in \"Cluster Data\" dialog", () => check(page, el("\"Show Scatter Plot\" input in \"Cluster Data\" dialog")), undefined, "Show Scatter Plot: see the clusters as you change the settings");
    await session.step(33, "Then the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
    await session.step(35, "When user clicks on OK button in \"Cluster Data\" dialog", () => clickOn(page, el("OK button in \"Cluster Data\" dialog")), undefined, "OK keeps them as a column");
    await session.step(36, "Then \"Cluster Data\" dialog should be hidden", () => shouldBe(page, el("\"Cluster Data\" dialog"), "hidden"));
    await session.step(37, "And the table should have a column \"clusters\"", () => hasColumn(page, "clusters"));
    await session.step(38, "And \"clusters\" column should have at least 2 distinct values", () => distinctValues(page, "clusters", 2));
    await session.step(40, "When user picks \"Petal.Width\" in the \"x\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "Petal.Width", "x", el("scatter plot viewer")), undefined, "By petal size, setosa forms a group of its own");
    await session.step(41, "Then \"x Column Name\" property of scatter plot viewer should be \"Petal.Width\"", () => propertyShouldBe(page, "x Column Name", el("scatter plot viewer"), "Petal.Width"));
    await session.step(42, "And no errors should have been logged", () => noErrors(page));
  });
});
