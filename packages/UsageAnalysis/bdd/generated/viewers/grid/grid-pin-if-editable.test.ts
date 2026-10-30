/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/grid/grid-pin-if-editable.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clearField, isExpanded, pressKeyIn, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {contextPanelOpen, contextPanelShows, enterOwnLogin, enterSharingLogin, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridPins} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Pin if editable pins a column for the users who may edit it", () => {
  const session = feature(test, "features/viewers/grid/grid-pin-if-editable.feature", import.meta.url);
  test("The switch pins the column for its editor, not for anyone else, and unpins it when off", {tag: ["@viewers", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(15, "And the grid should pin the columns \"\"", () => gridPins(page, ""));
    await session.step(16, "When user clicks on the \"header SEVERITY\" area of grid", () => clickArea(page, "header SEVERITY", el("grid")));
    await session.step(17, "Then the context panel should show \"SEVERITY\"", () => contextPanelShows(page, "SEVERITY"));
    await session.step(18, "Given Advanced accordion header in context panel is expanded", () => isExpanded(page, el("Advanced accordion header in context panel")));
    await session.step(19, "And Permissions accordion header in context panel is expanded", () => isExpanded(page, el("Permissions accordion header in context panel")));
    await session.step(22, "When user enters the current user's login into \"Editable by\" input in context panel", () => enterOwnLogin(page, el("\"Editable by\" input in context panel")));
    await session.step(23, "And user checks \"Pin if editable\" input in context panel", () => check(page, el("\"Pin if editable\" input in context panel")));
    await session.step(24, "Then the grid should pin the columns \"SEVERITY\"", () => gridPins(page, "SEVERITY"));
    await session.step(25, "When user enters the sharing user's login into \"Editable by\" input in context panel", () => enterSharingLogin(page, el("\"Editable by\" input in context panel")));
    await session.step(26, "Then the grid should pin the columns \"\"", () => gridPins(page, ""));
    await session.step(27, "When user enters the current user's login into \"Editable by\" input in context panel", () => enterOwnLogin(page, el("\"Editable by\" input in context panel")));
    await session.step(28, "Then the grid should pin the columns \"SEVERITY\"", () => gridPins(page, "SEVERITY"));
    await session.step(29, "When user unchecks \"Pin if editable\" input in context panel", () => uncheck(page, el("\"Pin if editable\" input in context panel")));
    await session.step(30, "Then the grid should pin the columns \"\"", () => gridPins(page, ""));
    await session.step(31, "When user clears \"Editable by\" input in context panel", () => clearField(page, el("\"Editable by\" input in context panel")));
    await session.step(32, "And user presses Enter in \"Editable by\" input in context panel", () => pressKeyIn(page, "Enter", el("\"Editable by\" input in context panel")));
    await session.step(33, "Then no errors should have been logged", () => noErrors(page));
  });
});
