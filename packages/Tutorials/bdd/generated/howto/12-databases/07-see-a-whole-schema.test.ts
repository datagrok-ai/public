/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/12-databases/07-see-a-whole-schema.feature
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
import {clickOn, expand, hoverOver, shouldBe, videoSays, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, contextPanelOpen, contextPanelShows, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Databases", () => {
  const session = feature(test, "features/howto/12-databases/07-see-a-whole-schema.feature", import.meta.url);
  test("See a whole schema at once", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And simple mode is off", () => simpleModeOff(page));
    await session.step(14, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(15, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(16, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(17, "When the video says \"One page with every table of a schema\"", () => videoSays(page, "One page with every table of a schema"));
    await session.step(18, "And user expands Databases tree node inside browse tree", () => expand(page, el("Databases tree node inside browse tree")));
    await session.step(19, "And user expands Databases---PostgresDart tree node inside browse tree", () => expand(page, el("Databases---PostgresDart tree node inside browse tree")));
    await session.step(20, "And user expands Databases---PostgresDart---Biologics tree node inside browse tree", () => expand(page, el("Databases---PostgresDart---Biologics tree node inside browse tree")));
    await session.step(21, "And user expands Databases---PostgresDart---Biologics---Schemas tree node inside browse tree", () => expand(page, el("Databases---PostgresDart---Biologics---Schemas tree node inside browse tree")));
    await session.step(23, "And user picks \"Browse\" from the context menu of Databases---PostgresDart---Biologics---Schemas---biologics tree node inside browse tree", () => pickFromContextMenu(page, "Browse", el("Databases---PostgresDart---Biologics---Schemas---biologics tree node inside browse tree")), undefined, "Browse on a schema opens the schema view");
    await session.step(24, "Then the \"Schema: biologics\" view should be current", () => viewIsCurrent(page, "Schema: biologics"));
    await session.step(25, "And \"drugs\" schema table should be visible", () => shouldBe(page, el("\"drugs\" schema table"), "visible"));
    await session.step(27, "When user hovers over \"expression_batches\" schema table", () => hoverOver(page, el("\"expression_batches\" schema table")), undefined, "Each box is a table and its columns");
    await session.step(29, "And user clicks on \"sequences\" text inside \"sequences\" schema table", () => clickOn(page, el("\"sequences\" text inside \"sequences\" schema table")), undefined, "Click a table's name for its details");
    await session.step(30, "Then the context panel should show \"sequences\"", () => contextPanelShows(page, "sequences"));
    await session.step(32, "When user clicks on \"heavy_chain\" text inside \"sequences\" schema table", () => clickOn(page, el("\"heavy_chain\" text inside \"sequences\" schema table")), undefined, "Click a column for the column's details");
    await session.step(33, "Then the context panel should show \"heavy_chain\"", () => contextPanelShows(page, "heavy_chain"));
    await session.step(35, "When user picks \"Get Top 100\" from the context menu of \"drugs\" schema table", () => pickFromContextMenu(page, "Get Top 100", el("\"drugs\" schema table")), undefined, "Right-click a box for the table's commands");
    await session.step(36, "Then the \"drugs\" view should be current", () => viewIsCurrent(page, "drugs"));
    await session.step(37, "And the table should have 100 rows", () => rowCount(page, 100));
    await session.step(38, "When the video says \"Get Top 100: the first hundred rows of the table\" at the bottom", () => videoSaysAt(page, "Get Top 100: the first hundred rows of the table", "bottom"));
    await session.step(39, "And user hovers over \"Rows: 100\" text in status bar", () => hoverOver(page, el("\"Rows: 100\" text in status bar")));
    await session.step(40, "Then no errors should have been logged", () => noErrors(page));
  });
});
