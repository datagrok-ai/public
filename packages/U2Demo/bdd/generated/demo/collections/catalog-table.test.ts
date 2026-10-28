/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demo/collections/catalog-table.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.catalog-table]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {openDemoPage} from '../../../bindings/demo.js';
import {catalogChoicesCleared, hoverAsGrid} from '../../../bindings/steps.js';
import {clickOn, collapse, doubleClickOn, expand, hoverOver, rightClickOn, shouldBe, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {itemShouldBe, shouldHaveItemsInTotal, shouldListInOrder, shouldReportReading} from '@datagrok-libraries/bdd/bindings/common/u2-steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Catalog table", () => {
  const session = feature(test, "features/demo/collections/catalog-table.feature", import.meta.url);
  test("A double click on a header sorts up, down, then back to the default order", {tag: ["@demo", "@realizes:u2.catalog-table"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given the catalog table choices kept under \"u2demo-campaigns\" are cleared now and put back at feature end", () => catalogChoicesCleared(page, "u2demo-campaigns"));
    await session.step(10, "And user opens the \"Catalog table\" demo page", () => openDemoPage(page, "Catalog table"));
    enter(page, "U2 Demo");
    await session.step(13, "Then campaigns catalog table should have 50 items in total", () => shouldHaveItemsInTotal(page, el("campaigns catalog table"), 50));
    await session.step(14, "And campaigns catalog table should report \"sort\" as \"\"", () => shouldReportReading(page, el("campaigns catalog table"), "sort", ""));
    await session.step(15, "When user double-clicks \"Molecules\" column header in campaigns catalog table", () => doubleClickOn(page, el("\"Molecules\" column header in campaigns catalog table")));
    await session.step(16, "Then campaigns catalog table should report \"sort\" as \"Molecules asc\"", () => shouldReportReading(page, el("campaigns catalog table"), "sort", "Molecules asc"));
    await session.step(17, "And campaigns catalog table should list the items \"KHL-4, KHL-14, SKP-9\" in this order", () => shouldListInOrder(page, el("campaigns catalog table"), "KHL-4, KHL-14, SKP-9"));
    await session.step(18, "When user double-clicks \"Molecules\" column header in campaigns catalog table", () => doubleClickOn(page, el("\"Molecules\" column header in campaigns catalog table")));
    await session.step(19, "Then campaigns catalog table should report \"sort\" as \"Molecules desc\"", () => shouldReportReading(page, el("campaigns catalog table"), "sort", "Molecules desc"));
    await session.step(20, "And campaigns catalog table should list the items \"SKP-9, KHL-14, KHL-4\" in this order", () => shouldListInOrder(page, el("campaigns catalog table"), "SKP-9, KHL-14, KHL-4"));
    await session.step(21, "And value of sort readout should have text \"Molecules descending\"", () => shouldHaveText(page, el("value of sort readout"), "Molecules descending"));
    await session.step(22, "When user double-clicks \"Molecules\" column header in campaigns catalog table", () => doubleClickOn(page, el("\"Molecules\" column header in campaigns catalog table")));
    await session.step(23, "Then campaigns catalog table should report \"sort\" as \"\"", () => shouldReportReading(page, el("campaigns catalog table"), "sort", ""));
    await session.step(24, "And value of sort readout should have text \"(last modified)\"", () => shouldHaveText(page, el("value of sort readout"), "(last modified)"));
  });
  test("The group menu groups the rows into counted sections that fold", {tag: ["@demo", "@realizes:u2.catalog-table"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given the catalog table choices kept under \"u2demo-campaigns\" are cleared now and put back at feature end", () => catalogChoicesCleared(page, "u2demo-campaigns"));
    await session.step(10, "And user opens the \"Catalog table\" demo page", () => openDemoPage(page, "Catalog table"));
    enter(page, "U2 Demo");
    await session.step(27, "When user clicks on campaigns-group icon button", () => clickOn(page, el("campaigns-group icon button")));
    await session.step(28, "And user clicks on \"Template\" menu item", () => clickOn(page, el("\"Template\" menu item")));
    await session.step(29, "Then campaigns catalog table should report \"groupBy\" as \"Template\"", () => shouldReportReading(page, el("campaigns catalog table"), "groupBy", "Template"));
    await session.step(30, "And campaigns catalog table should report \"groups\" as \"Fragment expansion (10), Quick sketchpad (17), Kinase hit-to-lead (14), GPCR series (9)\"", () => shouldReportReading(page, el("campaigns catalog table"), "groups", "Fragment expansion (10), Quick sketchpad (17), Kinase hit-to-lead (14), GPCR series (9)"));
    await session.step(31, "And item \"Quick sketchpad\" of campaigns catalog table should be expanded", () => itemShouldBe(page, "Quick sketchpad", el("campaigns catalog table"), "expanded"));
    await session.step(32, "When user collapses \"Quick sketchpad\" table row in campaigns catalog table", () => collapse(page, el("\"Quick sketchpad\" table row in campaigns catalog table")));
    await session.step(33, "Then item \"Quick sketchpad\" of campaigns catalog table should be collapsed", () => itemShouldBe(page, "Quick sketchpad", el("campaigns catalog table"), "collapsed"));
    await session.step(34, "And campaigns catalog table should report \"collapsed\" as \"SKP\"", () => shouldReportReading(page, el("campaigns catalog table"), "collapsed", "SKP"));
    await session.step(35, "And campaigns catalog table should report \"rows\" as \"33\"", () => shouldReportReading(page, el("campaigns catalog table"), "rows", "33"));
    await session.step(36, "When user expands \"Quick sketchpad\" table row in campaigns catalog table", () => expand(page, el("\"Quick sketchpad\" table row in campaigns catalog table")));
    await session.step(37, "Then campaigns catalog table should report \"rows\" as \"50\"", () => shouldReportReading(page, el("campaigns catalog table"), "rows", "50"));
  });
  test("The row under the pointer is marked as the platform grid marks it, lighter than a section header", {tag: ["@demo", "@realizes:u2.catalog-table"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given the catalog table choices kept under \"u2demo-campaigns\" are cleared now and put back at feature end", () => catalogChoicesCleared(page, "u2demo-campaigns"));
    await session.step(10, "And user opens the \"Catalog table\" demo page", () => openDemoPage(page, "Catalog table"));
    enter(page, "U2 Demo");
    await session.step(40, "When user clicks on campaigns-group icon button", () => clickOn(page, el("campaigns-group icon button")));
    await session.step(41, "And user clicks on \"Template\" menu item", () => clickOn(page, el("\"Template\" menu item")));
    await session.step(42, "And user hovers over \"SKP-9\" table row in campaigns catalog table", () => hoverOver(page, el("\"SKP-9\" table row in campaigns catalog table")));
    await session.step(43, "Then the row under the pointer in campaigns catalog table should be marked as the platform grid marks it", () => hoverAsGrid(page, el("campaigns catalog table")));
  });
  test("The eye's menu adds a column at the end", {tag: ["@demo", "@realizes:u2.catalog-table"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given the catalog table choices kept under \"u2demo-campaigns\" are cleared now and put back at feature end", () => catalogChoicesCleared(page, "u2demo-campaigns"));
    await session.step(10, "And user opens the \"Catalog table\" demo page", () => openDemoPage(page, "Catalog table"));
    enter(page, "U2 Demo");
    await session.step(46, "When user clicks on campaigns-columns icon button", () => clickOn(page, el("campaigns-columns icon button")));
    await session.step(47, "And user clicks on \"Last Modified\" menu item", () => clickOn(page, el("\"Last Modified\" menu item")));
    await session.step(48, "Then campaigns catalog table should report \"columns\" as \"Name,Code,Created,Author,Last Modified by,Molecules,Status,Last Modified\"", () => shouldReportReading(page, el("campaigns catalog table"), "columns", "Name,Code,Created,Author,Last Modified by,Molecules,Status,Last Modified"));
    await session.step(49, "And \"Last Modified\" column header in campaigns catalog table should be visible", () => shouldBe(page, el("\"Last Modified\" column header in campaigns catalog table"), "visible"));
  });
  test("A row's actions show once its access is known, and are its context menu; a view-only row has none", {tag: ["@demo", "@realizes:u2.catalog-table"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given the catalog table choices kept under \"u2demo-campaigns\" are cleared now and put back at feature end", () => catalogChoicesCleared(page, "u2demo-campaigns"));
    await session.step(10, "And user opens the \"Catalog table\" demo page", () => openDemoPage(page, "Catalog table"));
    enter(page, "U2 Demo");
    await session.step(52, "Then \"Delete campaign\" icon button in \"SKP-9\" table row in campaigns catalog table should be visible", () => shouldBe(page, el("\"Delete campaign\" icon button in \"SKP-9\" table row in campaigns catalog table"), "visible"));
    await session.step(53, "And \"Delete campaign\" icon button in \"FRAG-9\" table row in campaigns catalog table should be absent", () => shouldBe(page, el("\"Delete campaign\" icon button in \"FRAG-9\" table row in campaigns catalog table"), "absent"));
    await session.step(54, "When user right-clicks on \"SKP-9\" table row in campaigns catalog table", () => rightClickOn(page, el("\"SKP-9\" table row in campaigns catalog table")));
    await session.step(55, "Then \"Delete campaign\" menu item should be visible", () => shouldBe(page, el("\"Delete campaign\" menu item"), "visible"));
    await session.step(56, "And \"Manage campaign permissions\" menu item should be visible", () => shouldBe(page, el("\"Manage campaign permissions\" menu item"), "visible"));
  });
});
