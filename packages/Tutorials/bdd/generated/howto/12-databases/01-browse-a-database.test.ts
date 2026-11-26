/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/12-databases/01-browse-a-database.feature
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
import {columnsExactly} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelOpen, contextPanelShows, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Databases", () => {
  const session = feature(test, "features/howto/12-databases/01-browse-a-database.feature", import.meta.url);
  test("Browse a database", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(16, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(17, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(18, "When the video says \"Databases are in Browse, under Databases\"", () => videoSays(page, "Databases are in Browse, under Databases"));
    await session.step(19, "And user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(20, "And user expands Databases tree node inside browse tree", () => expand(page, el("Databases tree node inside browse tree")));
    await session.step(22, "And user expands Databases---PostgresDart tree node inside browse tree", () => expand(page, el("Databases---PostgresDart tree node inside browse tree")), undefined, "Connections are grouped by database type");
    await session.step(24, "And user hovers over Databases---PostgresDart---Biologics tree node inside browse tree", () => hoverOver(page, el("Databases---PostgresDart---Biologics tree node inside browse tree")), undefined, "Point at a connection to see its server and database");
    await session.step(26, "And user expands Databases---PostgresDart---Biologics tree node inside browse tree", () => expand(page, el("Databases---PostgresDart---Biologics tree node inside browse tree")), undefined, "A connection holds its saved queries and its schemas");
    await session.step(27, "And user expands Databases---PostgresDart---Biologics---Schemas tree node inside browse tree", () => expand(page, el("Databases---PostgresDart---Biologics---Schemas tree node inside browse tree")));
    await session.step(29, "And user expands Databases---PostgresDart---Biologics---Schemas---biologics tree node inside browse tree", () => expand(page, el("Databases---PostgresDart---Biologics---Schemas---biologics tree node inside browse tree")), undefined, "A schema lists its tables");
    await session.step(31, "And user clicks on Databases---PostgresDart---Biologics---Schemas---biologics---drugs tree node inside browse tree", () => clickOn(page, el("Databases---PostgresDart---Biologics---Schemas---biologics---drugs tree node inside browse tree")), undefined, "Click a table to preview its first rows");
    await session.step(32, "Then the context panel should show \"drugs\"", () => contextPanelShows(page, "drugs"));
    await session.step(34, "When user expands Columns section in context panel", () => expand(page, el("Columns section in context panel")), undefined, "The context panel lists its columns");
    await session.step(36, "And user expands Databases---PostgresDart---Biologics---Schemas---biologics---drugs tree node inside browse tree", () => expand(page, el("Databases---PostgresDart---Biologics---Schemas---biologics---drugs tree node inside browse tree")), undefined, "Or expand the table in the tree");
    await session.step(37, "Then Databases---PostgresDart---Biologics---Schemas---biologics---drugs---smiles tree node inside browse tree should be visible", () => shouldBe(page, el("Databases---PostgresDart---Biologics---Schemas---biologics---drugs---smiles tree node inside browse tree"), "visible"));
    await session.step(39, "When user picks \"Get All\" from the context menu of Databases---PostgresDart---Biologics---Schemas---biologics---drugs tree node inside browse tree", () => pickFromContextMenu(page, "Get All", el("Databases---PostgresDart---Biologics---Schemas---biologics---drugs tree node inside browse tree")), undefined, "Get All opens the whole table");
    await session.step(40, "Then the \"drugs\" view should be current", () => viewIsCurrent(page, "drugs"));
    await session.step(41, "And the table should have 200 rows", () => rowCount(page, 200));
    await session.step(42, "And the table should have the columns \"id, identifier, name, smiles\"", () => columnsExactly(page, "id, identifier, name, smiles"));
    await session.step(43, "When the video says \"Now it is a table like any other: sort, filter, chart\" at the bottom", () => videoSaysAt(page, "Now it is a table like any other: sort, filter, chart", "bottom"));
    await session.step(44, "And user hovers over \"Rows: 200\" text in status bar", () => hoverOver(page, el("\"Rows: 200\" text in status bar")));
    await session.step(45, "Then no errors should have been logged", () => noErrors(page));
  });
});
