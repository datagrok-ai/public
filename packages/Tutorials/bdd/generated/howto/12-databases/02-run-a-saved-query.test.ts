/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/12-databases/02-run-a-saved-query.feature
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
import {clickOn, expand, hoverOver, selectIn, shouldBe, videoSays, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {everyValueContains} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {clickRunsQuery} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {browsePanelOpen, contextPanelClosed, currentViewType, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickFromContextMenu, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Databases", () => {
  const session = feature(test, "features/howto/12-databases/02-run-a-saved-query.feature", import.meta.url);
  test("Run a saved query", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(17, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(18, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(19, "When the video says \"Saved queries are listed under their connection\"", () => videoSays(page, "Saved queries are listed under their connection"));
    await session.step(20, "And user expands Databases tree node inside browse tree", () => expand(page, el("Databases tree node inside browse tree")));
    await session.step(21, "And user expands Databases---PostgresDart tree node inside browse tree", () => expand(page, el("Databases---PostgresDart tree node inside browse tree")));
    await session.step(22, "And user expands Databases---PostgresDart---Biologics tree node inside browse tree", () => expand(page, el("Databases---PostgresDart---Biologics tree node inside browse tree")));
    await session.step(24, "And user picks \"Run\" from the context menu of Databases---PostgresDart---Biologics---Antibody-Profile-by-Organism-and-Target tree node inside browse tree", () => pickFromContextMenu(page, "Run", el("Databases---PostgresDart---Biologics---Antibody-Profile-by-Organism-and-Target tree node inside browse tree")), undefined, "Run asks for the query's parameters");
    await session.step(25, "Then \"Antibody Profile by Organism and Target\" dialog should be visible", () => shouldBe(page, el("\"Antibody Profile by Organism and Target\" dialog"), "visible"));
    await session.step(27, "When user selects \"E. coli\" in Organism input in \"Antibody Profile by Organism and Target\" dialog", () => selectIn(page, "E. coli", el("Organism input in \"Antibody Profile by Organism and Target\" dialog")), undefined, "The choices come from the database");
    await session.step(28, "And user selects \"HER2\" in Target input in \"Antibody Profile by Organism and Target\" dialog", () => selectIn(page, "HER2", el("Target input in \"Antibody Profile by Organism and Target\" dialog")));
    await session.step(29, "And user clicks on OK button in \"Antibody Profile by Organism and Target\" dialog", () => clickOn(page, el("OK button in \"Antibody Profile by Organism and Target\" dialog")));
    await session.step(30, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
    await session.step(31, "And every value of \"organism_name\" column should contain \"E. coli\"", () => everyValueContains(page, "organism_name", "E. coli"));
    await session.step(32, "And every value of \"target_name\" column should contain \"HER2\"", () => everyValueContains(page, "target_name", "HER2"));
    await session.step(33, "And the open tableview should have 1 pie chart viewer", () => viewerCount(page, 1, "pie chart"));
    await session.step(34, "When the video says \"The result opens with the layout saved with the query\" at the top", () => videoSaysAt(page, "The result opens with the layout saved with the query", "top"));
    await session.step(35, "And user hovers over pie chart viewer", () => hoverOver(page, el("pie chart viewer")));
    await session.step(36, "Then toolbox should be visible", () => shouldBe(page, el("toolbox"), "visible"));
    await session.step(38, "When user selects \"PD-1\" in Target input in toolbox", () => selectIn(page, "PD-1", el("Target input in toolbox")), undefined, "The toolbox keeps the parameters: pick another target");
    await session.step(40, "And user clicks on REFRESH button in toolbox and the query it runs completes", () => clickRunsQuery(page, el("REFRESH button in toolbox")), undefined, "REFRESH runs the query again");
    await session.step(41, "Then every value of \"target_name\" column should contain \"PD-1\"", () => everyValueContains(page, "target_name", "PD-1"));
    await session.step(42, "And every value of \"organism_name\" column should contain \"E. coli\"", () => everyValueContains(page, "organism_name", "E. coli"));
    await session.step(43, "When the video says \"Same view, new rows from the database\" at the top", () => videoSaysAt(page, "Same view, new rows from the database", "top"));
    await session.step(44, "And user hovers over pie chart viewer", () => hoverOver(page, el("pie chart viewer")));
    await session.step(45, "Then no errors should have been logged", () => noErrors(page));
  });
});
