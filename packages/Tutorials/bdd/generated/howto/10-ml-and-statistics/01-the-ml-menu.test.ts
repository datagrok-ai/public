/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/10-ml-and-statistics/01-the-ml-menu.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
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
import {pointAtTopMenuItem} from '../../../bindings/howto-data.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · ML", () => {
  const session = feature(test, "features/howto/10-ml-and-statistics/01-the-ml-menu.feature", import.meta.url);
  test("What the ML menu does", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(19, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(20, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(21, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(22, "When the video says \"ML: statistics and machine learning on the open table\" at the top", () => videoSaysAt(page, "ML: statistics and machine learning on the open table", "top"));
    await session.step(24, "And user points at \"ML > Analyze > PCA...\" in the top menu", () => pointAtTopMenuItem(page, "ML > Analyze > PCA..."), undefined, "Analyze: PCA, PLS and multivariate analysis");
    await session.step(25, "And user points at \"ML > Analyze > Multivariate Analysis...\" in the top menu", () => pointAtTopMenuItem(page, "ML > Analyze > Multivariate Analysis..."));
    await session.step(27, "And user points at \"ML > Analyze > Group Comparison\" in the top menu", () => pointAtTopMenuItem(page, "ML > Analyze > Group Comparison"), undefined, "Group Comparison: t-test, ANOVA, each group against a control");
    await session.step(28, "And user points at \"ML > Analyze > Group Comparison > ANOVA...\" in the top menu", () => pointAtTopMenuItem(page, "ML > Analyze > Group Comparison > ANOVA..."));
    await session.step(30, "And user points at \"ML > Reduce Dimensionality...\" in the top menu", () => pointAtTopMenuItem(page, "ML > Reduce Dimensionality..."), undefined, "Reduce Dimensionality: a 2D map of the rows (UMAP, t-SNE)");
    await session.step(32, "And user points at \"ML > Impute Missing Values...\" in the top menu", () => pointAtTopMenuItem(page, "ML > Impute Missing Values..."), undefined, "Impute Missing Values: fill the gaps from similar rows");
    await session.step(34, "And user points at \"ML > Pareto Front...\" in the top menu", () => pointAtTopMenuItem(page, "ML > Pareto Front..."), undefined, "Pareto Front: the best trade-offs between several goals");
    await session.step(36, "And user points at \"ML > Cluster\" in the top menu", () => pointAtTopMenuItem(page, "ML > Cluster"), undefined, "Cluster: k-means, DBSCAN, MCL and a hierarchical tree");
    await session.step(37, "And user points at \"ML > Cluster > Cluster Data...\" in the top menu", () => pointAtTopMenuItem(page, "ML > Cluster > Cluster Data..."));
    await session.step(39, "And user points at \"ML > Models\" in the top menu", () => pointAtTopMenuItem(page, "ML > Models"), undefined, "Models: train a model on this table, apply it to another");
    await session.step(40, "And user points at \"ML > Models > Train Model...\" in the top menu", () => pointAtTopMenuItem(page, "ML > Models > Train Model..."));
    await session.step(41, "Then \"Train Model...\" menu item should be visible", () => shouldBe(page, el("\"Train Model...\" menu item"), "visible"));
    await session.step(42, "And no errors should have been logged", () => noErrors(page));
  });
});
