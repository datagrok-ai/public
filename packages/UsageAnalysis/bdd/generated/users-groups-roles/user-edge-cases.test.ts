/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/users-groups-roles/user-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.users]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, expand, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, galleryCountLower, rememberGalleryCount, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, menuDoesNotList, menuLists, noBalloons, noErrors, openContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("User edge cases: the running account cannot disable itself from the Users view", () => {
  const session = feature(test, "features/users-groups-roles/user-edge-cases.feature", import.meta.url);
  test("The running account's own context menu offers no Disable...", {tag: ["@users", "@edge-cases", "@realizes:views.users"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(15, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await session.step(16, "And user clicks on \"Platform > Users\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Users\" tree node inside browse tree")));
    await session.step(19, "Then the \"Users\" view should be current", () => viewIsCurrent(page, "Users"));
    await session.step(20, "When user remembers the gallery counter", () => rememberGalleryCount(page));
    await session.step(21, "And user types \"{login}\" into gallery search", () => typeInto(page, session.text("{login}"), el("gallery search")));
    await session.step(22, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
    await session.step(23, "When user opens the context menu of \"{login}\" link in gallery", () => openContextMenu(page, el(session.text("\"{login}\" link in gallery"))));
    await session.step(24, "Then the open menu should list \"Details\"", () => menuLists(page, "Details"));
    await session.step(25, "And the open menu should not list \"Disable...\"", () => menuDoesNotList(page, "Disable..."));
    await session.step(26, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(27, "And user clears gallery search", () => clearField(page, el("gallery search")));
    await session.step(28, "Then no errors should have been logged", () => noErrors(page));
    await session.step(29, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
