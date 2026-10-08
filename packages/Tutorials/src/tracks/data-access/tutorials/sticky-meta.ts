import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import {Tutorial, TutorialPrerequisites} from '@datagrok-libraries/tutorials/src/tutorial';
import {getPlatform, Platform} from '../../shortcuts';
import {_package} from '../../../package';
import {fromEvent, interval} from 'rxjs';
import {filter} from 'rxjs/operators';
import {describeElements} from '../../compute/tutorials/utils';
import {elementClick} from '../../eda/tutorials/utils';

enum LINKS {
  STICKY_META = 'https://datagrok.ai/help/govern/catalog/sticky-meta',
}

/** What the steps tell the learner to type. The later steps look these values back up, so they have
 * to be the same string in both places — they used to disagree in case, and the lookup that missed
 * took the rest of the tutorial down with it. */
const SCHEMA_NAME = 'schema for tutorial';
const PROPERTY_NAME = 'project name';
const PROPERTY_VALUE = 'Tutorial';
const TABLE_NAME = 'Sticky Meta molecules';

/** [Tutorial.waitFor], naming the element that was missing when it never turns up — a lookup that
 * silently returns `undefined` throws a TypeError one line later. */
async function waitForElement<T extends Element>(get: () => T | null | undefined, what: string): Promise<T> {
  const el = await Tutorial.waitFor(get);
  if (el == null)
    throw new Error(`Sticky Meta tutorial: ${what} never appeared`);
  return el;
}

/** Case-insensitive match on an element's trimmed text. */
const textIs = (value: string) => (el: Element) => el.textContent?.trim().toLowerCase() === value.toLowerCase();

/** Sticky Meta tutorial */
export class StickyMetaTutorial extends Tutorial {
  get name() {return 'Sticky Meta';}
  get description(): string {
    return `Learn how to create and use Sticky Meta in Datagrok. 
      Define schemas and entity types, annotate data, and explore metadata across datasets.`;
  }
  get steps() {return 20;}
  get icon() {return '📌';}

  helpUrl: string = LINKS.STICKY_META;
  prerequisites: TutorialPrerequisites = {packages: ['Chem']};
  demoTable: string = '';
  platform: Platform = getPlatform();

  protected async _run() {
    grok.shell.windows.showContextPanel = false;
    grok.shell.windows.showProjects = false;
    grok.shell.windows.showToolbox = false;
    grok.shell.windows.showBrowse = true;

    this.header.textContent = this.name;

    // --- Step 1: Introduction ---
    this.title('Introduction', true);
    this.describe('Sticky Meta lets you annotate data with structured metadata. ' +
      'This metadata is stored centrally and is integrated across Datagrok for search, filtering, and analysis.');
    this.describe(ui.link('Learn more', this.helpUrl).outerHTML);

    // --- Step 2: Open types node ---
    this.title('Explore entity types');
    // the Sticky Meta group loads its children after it expands, so the nodes are looked up on every poll
    const stickyMetaNode = () => {
      const platformNode = findTreeNode(grok.shell.browsePanel.mainTree, 'Platform', true);
      return platformNode ? findTreeNode(platformNode, 'Sticky Meta', true) : undefined;
    };
    const stickyMetaChild = (name: string) => {
      const group = stickyMetaNode();
      // the row, not its caption: the platform opens a node on a click anywhere in it
      return group ? findTreeNode(group, name, true)?.root ?? null : null;
    };
    const button = (caption: string) => Array.from(document.querySelectorAll('button'))
      .find((b) => b.textContent?.trim() === caption) ?? null;
    stickyMetaNode();

    this.describe('Entity types define the objects you annotate, e.g., molecules. ' +
      'They are located under <b>Browse → Platform → Sticky Meta → Types.</b>');
    await this.action(
      'Open Types node',
      elementClick(() => stickyMetaChild('Types')),
      () => stickyMetaChild('Types'),
      '',
      () => stickyMetaChild('Types')!.click(),
    );

    // --- Step 3: Open new entity type dialog ---
    await Tutorial.waitFor(() => button('New Entity Type...'));
    const typeDialog = await this.openDialog(
      'Create a new entity type',
      'Create a new entity type',
      () => button('New Entity Type...'),
      'Click "New Entity Type..." to create a new entity type.',
      () => button('New Entity Type...')!.click(),
    );

    // --- Step 4: Explore entity type dialog ---
    const typeDialogRoot = typeDialog.root;
    const nameInput = Array.from(typeDialogRoot.querySelectorAll('label.ui-label span'))
      .find((el) => el.textContent?.trim() === 'Name') as HTMLElement;
    const matchInput = Array.from(typeDialogRoot.querySelectorAll('label.ui-label span'))
      .find((el) => el.textContent?.trim() === 'Matching expression') as HTMLElement;

    let doneBtn = describeElements([nameInput, matchInput], [
      '# Name\nEntity type name, e.g., "Molecule".',
      '# Matching expression\nDefines objects in this type, e.g., "semtype=molecule".',
    ]);

    await this.action(
      'Explore entity type dialog',
      fromEvent(doneBtn, 'click'),
      undefined,
      'Click "Next" to proceed.',
      () => doneBtn.click(),
    );

    // --- Step 5: Fill in entity type fields ---
    await this.dlgInputAction(typeDialog, 'Set "Name" to "molecule-tutorial"', 'Name', 'molecule-tutorial');
    await this.dlgInputAction(typeDialog, 'Set "Matching expression" to "semtype=Molecule"', 'Matching expression', 'semtype=Molecule');

    await this.dialogOkAction(typeDialog, 'Save entity type', 'Click OK to save the entity type.');

    // --- Step 6: Open schemas node ---
    this.title('Explore schemas');
    this.describe('Schemas define the metadata fields and are linked to entity types. ' +
      'They are located under <b>Browse → Platform → Sticky Meta → Schemas.</b>');
    await this.action(
      'Open schemas node',
      elementClick(() => stickyMetaChild('Schemas')),
      () => stickyMetaChild('Schemas'),
      '',
      () => stickyMetaChild('Schemas')!.click(),
    );

    // --- Step 7: Open new schema dialog ---
    await Tutorial.waitFor(() => button('New Schema...'));
    const schemaDialog = await this.openDialog(
      'Create a new schema',
      'Create a new schema',
      () => button('New Schema...'),
      'Click "New Schema..." to create a schema.',
      () => button('New Schema...')!.click(),
    );

    // --- Step 8: Explore schema dialog ---
    const schemaDialogRoot = schemaDialog.root;
    const schemaNameInput = Array.from(schemaDialogRoot.querySelectorAll('.ui-input-root'))
      .find((el) => el.querySelector('label.ui-label span')?.textContent?.trim() === 'Name')
      ?.querySelector('input') as HTMLInputElement;
    const assocWithLabel = Array.from(schemaDialogRoot.querySelectorAll('.ui-input-root'))
      .find((el) => el.querySelector('label.ui-label span')?.textContent?.trim() === 'Associated with:')
      ?.querySelector('div.ui-input-editor') as HTMLDivElement;
    const propertyTable = schemaDialogRoot.querySelector('table.d4-item-table') as HTMLTableElement;
    const propertyRow = propertyTable?.querySelectorAll('tbody > tr')[1] as HTMLTableRowElement;

    doneBtn = describeElements([schemaNameInput, assocWithLabel, propertyRow], [
      '# Name\nName of the schema.',
      '# Associated with\nThe entity type this schema applies to.',
      '# Properties\nMetadata fields with a name and type (string, int, bool, double, datetime).',
    ]);

    await this.action('Explore schema dialog', fromEvent(doneBtn, 'click'), null, '', () => doneBtn.click());

    // --- Step 9: Fill in schema details ---
    await this.textInpAction(schemaDialogRoot, `Set "Name" to "${SCHEMA_NAME}"`, 'Name', SCHEMA_NAME);

    const assocWithRoot = $(schemaDialogRoot)
      .find('label.ui-label.ui-input-label span')
      .filter((_, el) => el.textContent?.trim() === 'Associated with:')
      .closest('.ui-input-root')[0]!;

    const selectEntitiesLabel = assocWithRoot.querySelector('.d4-link-action') as HTMLElement;
    await this.action(
      'Select associated entity',
      new Promise<void>((resolve) => selectEntitiesLabel.addEventListener('click', () => resolve(), {once: true})),
      selectEntitiesLabel,
      '',
      () => selectEntitiesLabel.click(),
    );

    const selectorDlg = await new Promise<DG.Dialog>((resolve) => {
      const iv = setInterval(() => {
        const dlg = DG.Dialog.getOpenDialogs().find((d) => d.title.includes('Select types for'));
        if (dlg) {clearInterval(iv); resolve(dlg);}
      }, 50);
    });

    // the rows are rendered after the dialog opens, so the checkbox is looked up on every poll
    const moleculeCheckbox = () => Array.from(selectorDlg.root.querySelectorAll('.property-grid-item-name-text'))
      .find((x) => x.textContent?.trim() === 'molecule-tutorial')
      ?.closest('.property-grid-item')?.querySelector('input[type="checkbox"]') as HTMLInputElement ?? null;

    await this.action(
      'Select molecule-tutorial',
      interval(200).pipe(filter(() => moleculeCheckbox()?.checked === true)),
      moleculeCheckbox,
      '',
      () => moleculeCheckbox()!.click(),
    );

    await this.dialogOkAction(selectorDlg, 'Confirm entity selection');

    await this.textInpAction(propertyRow, `Set property "Name" to "${PROPERTY_NAME}"`, 'Name', PROPERTY_NAME);
    await this.choiceInputAction(propertyRow, 'Set property "Type" to "string"', 'Property Type', 'string');
    await this.dialogOkAction(schemaDialog, 'Save schema', 'Click OK to save schema.');

    // --- Step 10: Open dataset ---
    this.title('Annotate dataset');
    this.t = await grok.data.loadTable(`${_package.webRoot}files/smiles.csv`);
    this.t.name = TABLE_NAME;
    const tv = grok.shell.addTableView(this.t);
    const grid = tv.grid;

    const rowIndex = 0;
    const colName = 'smiles';
    this.describe(
      `The <b>${TABLE_NAME}</b> table has just opened. Blue circles in its cells indicate metadata availability. ` +
      `In this table, click the <b>first cell</b> in the <b>${colName}</b> column, ` +
      'then annotate this molecule in the <b>Sticky meta</b> section of the Context Panel. ' +
      'Stay in this table: annotations made in other tables do not count for this tutorial.',
    );
    grok.shell.windows.showContextPanel = true;

    await this.action(
      `In the ${TABLE_NAME} table, click the first cell in the ${colName} column`,
      new Promise<void>((resolve) => {
        const sub = grid.onCellClick.subscribe((args) => {
          if (args.cell.rowIndex === rowIndex && args.cell.column.name === colName) {
            sub.unsubscribe();
            resolve();
          }
        });
      }),
      null,
      '',
      () => Tutorial.clickCell(grid, colName, rowIndex),
    );

    // --- Step 11: Fill in sticky meta property ---
    grok.shell.windows.showContextPanel = true;
    const stickyPane = await waitForElement<HTMLElement>(
      () => document.querySelector<HTMLElement>('.d4-accordion-pane[d4-title="Sticky meta"]'),
      'the Sticky meta pane');
    const header = stickyPane.querySelector('.d4-accordion-pane-header') as HTMLElement;
    if (header && !header.classList.contains('expanded')) header.click();

    const schemaSection = await waitForElement<HTMLElement>(
      () => Array.from(stickyPane.querySelectorAll<HTMLElement>('.grok-sticky-meta-schema-header'))
        .find(textIs(SCHEMA_NAME)), `the "${SCHEMA_NAME}" section`);
    schemaSection.scrollIntoView({block: 'start'});

    // by its label, not by a generated input name: the property is named by the learner
    const projectPropertyInput = await waitForElement<HTMLInputElement>(
      () => Array.from(schemaSection.parentElement!.querySelectorAll('.ui-input-root'))
        .find((root) => textIs(PROPERTY_NAME)(root.querySelector('label.ui-label span')!))
        ?.querySelector('input'), `the "${PROPERTY_NAME}" input`);
    await this.action(
      `Set "${PROPERTY_NAME}" to "${PROPERTY_VALUE}"`,
      new Promise<void>((resolve) => {
        const listener = () => {
          if (projectPropertyInput.value.trim().toLowerCase() === PROPERTY_VALUE.toLowerCase()) {
            projectPropertyInput.removeEventListener('input', listener);
            resolve();
          }
        };
        projectPropertyInput.addEventListener('input', listener);
      }),
      [schemaSection, projectPropertyInput],
      `In the Context Panel, under <b>Sticky meta → ${SCHEMA_NAME}</b>, ` +
        `type <b>${PROPERTY_VALUE}</b> into the <b>${PROPERTY_NAME}</b> field.`,
      () => Tutorial.setInputValue(projectPropertyInput, PROPERTY_VALUE),
    );

    // --- Step 12: Save sticky meta ---
    const schemaSaveButton = () => {
      for (let el = schemaSection.nextElementSibling; el != null; el = el.nextElementSibling) {
        const btn = el.matches('button[name^="button-Save"]') ? el : el.querySelector('button[name^="button-Save"]');
        if (btn != null)
          return btn as HTMLButtonElement;
      }
      return null;
    };
    const saveBtn = await waitForElement<HTMLButtonElement>(() => {
      const btn = schemaSaveButton();
      return btn != null && !btn.classList.contains('disabled') && btn.getAttribute('aria-disabled') !== 'true' ?
        btn : null;
    }, `the enabled Save button of "${SCHEMA_NAME}"`);

    await this.action(
      `Click SAVE under "${SCHEMA_NAME}"`,
      new Promise<void>((resolve) => saveBtn.addEventListener('click', () => resolve(), {once: true})),
      saveBtn,
      `Click the <b>SAVE</b> button right below the <b>${PROPERTY_NAME}</b> field. ` +
        'Each schema in the pane has its own SAVE button.',
      () => saveBtn.click(),
    );

    await this.action(
      'Hover a cell to verify metadata tooltip.',
      new Promise<void>((resolve) => {
        const sub = grid.onCellTooltip((cell) => {
          if (cell.cell.column.name === 'smiles' && cell.gridRow === 0) {
            sub.unsubscribe();
            resolve();
          }
        });
      }),
    );

    this.title('Sticky Meta tutorial completed', true);
    this.describe('You have created an entity type and schema, annotated a dataset, and explored Sticky Meta.');
  }
}

/** Helper function to find nodes in tree view */
function findTreeNode(parent: DG.TreeViewGroup, name: string, expand: boolean = false): DG.TreeViewGroup | undefined {
  const node = parent.children.find((child) => child.text === name) as DG.TreeViewGroup | undefined;
  if (node && expand) node.expanded = true;
  return node;
}
