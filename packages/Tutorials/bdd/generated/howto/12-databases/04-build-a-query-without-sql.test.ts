/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/12-databases/04-build-a-query-without-sql.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
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
import {addToVisualQueryRow, visualQueryPreviewRows, visualQueryRowHolds} from '../../../bindings/howto-databases.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {expand, hoverOver, videoSays, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelClosed, currentViewType, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Databases", () => {
  const session = feature(test, "features/howto/12-databases/04-build-a-query-without-sql.feature", import.meta.url);
  test("Build a query without SQL", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(17, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(18, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(19, "When the video says \"A query built by picking columns, no SQL\"", () => videoSays(page, "A query built by picking columns, no SQL"));
    await session.step(20, "And user expands Databases tree node inside browse tree", () => expand(page, el("Databases tree node inside browse tree")));
    await session.step(21, "And user expands Databases---PostgresDart tree node inside browse tree", () => expand(page, el("Databases---PostgresDart tree node inside browse tree")));
    await session.step(22, "And user expands Databases---PostgresDart---Biologics tree node inside browse tree", () => expand(page, el("Databases---PostgresDart---Biologics tree node inside browse tree")));
    await session.step(23, "And user expands Databases---PostgresDart---Biologics---Schemas tree node inside browse tree", () => expand(page, el("Databases---PostgresDart---Biologics---Schemas tree node inside browse tree")));
    await session.step(24, "And user expands Databases---PostgresDart---Biologics---Schemas---biologics tree node inside browse tree", () => expand(page, el("Databases---PostgresDart---Biologics---Schemas---biologics tree node inside browse tree")));
    await session.step(26, "And user picks \"New Visual Query...\" from the context menu of Databases---PostgresDart---Biologics---Schemas---biologics---expression_batches tree node inside browse tree", () => pickFromContextMenu(page, "New Visual Query...", el("Databases---PostgresDart---Biologics---Schemas---biologics---expression_batches tree node inside browse tree")), undefined, "New Visual Query... on a table opens the builder");
    await session.step(27, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
    await session.step(28, "When the video says \"The result under the builder follows every change\" at the bottom", () => videoSaysAt(page, "The result under the builder follows every change", "bottom"));
    await session.step(29, "And user hovers over \"Grid\" grid", () => hoverOver(page, el("\"Grid\" grid")));
    await session.step(31, "And user adds \"expression_system\" to the \"Group by\" row of the visual query", () => addToVisualQueryRow(page, "expression_system", "Group by"), undefined, "Group by: one row per expression system");
    await session.step(32, "Then the \"Group by\" row of the visual query should hold \"expression_system\"", () => visualQueryRowHolds(page, "Group by", "expression_system"));
    await session.step(33, "And the visual query should preview 3 rows", () => visualQueryPreviewRows(page, 3));
    await session.step(35, "When user adds \"yield_mg\" to the \"Aggregate\" row of the visual query", () => addToVisualQueryRow(page, "yield_mg", "Aggregate"), undefined, "Aggregate: the average yield of each");
    await session.step(36, "Then the \"Aggregate\" row of the visual query should hold \"avg(yield_mg)\"", () => visualQueryRowHolds(page, "Aggregate", "avg(yield_mg)"));
    await session.step(38, "When user adds \"chain\" to the \"Pivot\" row of the visual query", () => addToVisualQueryRow(page, "chain", "Pivot"), undefined, "Pivot: a column for each chain");
    await session.step(39, "Then the \"Pivot\" row of the visual query should hold \"chain\"", () => visualQueryRowHolds(page, "Pivot", "chain"));
    await session.step(40, "And the visual query should preview 3 rows", () => visualQueryPreviewRows(page, 3));
    await session.step(41, "When the video says \"Average yield by expression system and chain\" at the bottom", () => videoSaysAt(page, "Average yield by expression system and chain", "bottom"));
    await session.step(42, "And user hovers over \"Grid\" grid", () => hoverOver(page, el("\"Grid\" grid")));
    await session.step(43, "Then no errors should have been logged", () => noErrors(page));
  });
});
