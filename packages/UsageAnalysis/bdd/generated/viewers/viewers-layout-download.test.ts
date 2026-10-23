/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/viewers-layout-download.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.chrome]
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
import {downloadedContains, fileMatchingDownloaded, watchDownloads} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {openDataset, viewHoldsViewers} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dockToViewEdge, dockedAtViewEdge} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The layout of a docked view downloads as a .layout file", () => {
  const session = feature(test, "features/viewers/viewers-layout-download.feature", import.meta.url);
  test("View > Layout > Download saves the view's layout with its viewers", {tag: ["@viewers", "@realizes:viewers.chrome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a histogram viewer", () => addViewer(page, "histogram"));
    await session.step(14, "And user adds a bar chart viewer", () => addViewer(page, "bar chart"));
    await session.step(15, "Then the current view should hold at least 2 viewers", () => viewHoldsViewers(page, 2));
    await session.step(18, "When user docks histogram viewer to the right edge of the view", () => dockToViewEdge(page, el("histogram viewer"), "right"));
    await session.step(19, "Then histogram viewer should be docked along the right edge of the view", () => dockedAtViewEdge(page, el("histogram viewer"), "right"));
    await session.step(20, "Given user watches downloads", () => watchDownloads(page));
    await session.step(21, "When user picks \"View > Layout > Download\" from the top menu", () => pickFromTopMenu(page, "View > Layout > Download"));
    await session.step(22, "Then a file matching \"\\.layout$\" should have been downloaded", () => fileMatchingDownloaded(page, "\\.layout$"));
    await session.step(23, "And the downloaded file should contain \"Histogram\"", () => downloadedContains(page, "Histogram"));
    await session.step(24, "And the downloaded file should contain \"Bar chart\"", () => downloadedContains(page, "Bar chart"));
    await session.step(25, "And no errors should have been logged", () => noErrors(page));
  });
});
