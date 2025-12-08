import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import {Tutorial, TutorialPrerequisites} from '@datagrok-libraries/tutorials/src/tutorial';
import {getPlatform, Platform} from '../../shortcuts';
import {_package} from '../../../package';
import {fromEvent} from 'rxjs';
import {describeElements} from '../../compute/tutorials/utils';
import {waitForElementClick} from '../../eda/tutorials/utils';

enum LINKS {
  STICKY_META = 'https://datagrok.ai/help/govern/catalog/sticky-meta',
}

/** Sticky Meta tutorial */
export class StickyMetaTutorial extends Tutorial {
  get name() { return 'Sticky Meta'; }
  get description(): string {
    return `Learn how to create and use Sticky Meta in Datagrok. 
      Define schemas and entity types, annotate data, and explore metadata across datasets.`;
  }
  get steps() { return 19; }
  get icon() { return '📌'; }

  helpUrl: string = LINKS.STICKY_META;
  prerequisites: TutorialPrerequisites = { packages: ['Chem'] };
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
    const platformNode = findTreeNode(grok.shell.browsePanel.mainTree, 'Platform', true);
    const stickyMetaNode = platformNode ? findTreeNode(platformNode, 'Sticky Meta', true) : undefined;
    const typesNode = stickyMetaNode ? findTreeNode(stickyMetaNode, 'Types', true) : undefined;

    this.describe('Entity types define the objects you annotate, e.g., molecules. ' +
      'They are located under Browse → Platform → Sticky Meta → Types.');
    await this.action(
      'Open Types node',
      waitForElementClick(typesNode!.captionLabel),
      typesNode!.captionLabel
    );
    await new Promise(resolve => setTimeout(resolve, 500));

    // --- Step 3: Open new entity type dialog ---
    this.describe('Click "New Entity Type..." to create a new entity type.');
    const newEntityTypeButton = Array.from(document.querySelectorAll('button'))
      .find(b => b.textContent?.trim() === 'New Entity Type...');
    const typeDialog = await this.openDialog(
      'Create a new entity type',
      'Create a new entity type',
      newEntityTypeButton,
      'Click "New Entity Type..." to start'
    );

    // --- Step 4: Explore entity type dialog ---
    this.describe('Fill in the name and matching expression.');
    const typeDialogRoot = typeDialog.root;
    const nameInput = Array.from(typeDialogRoot.querySelectorAll('label.ui-label span'))
      .find(el => el.textContent?.trim() === 'Name') as HTMLElement;
    const matchInput = Array.from(typeDialogRoot.querySelectorAll('label.ui-label span'))
      .find(el => el.textContent?.trim() === 'Matching expression') as HTMLElement;

    let doneBtn = describeElements([nameInput, matchInput], [
      '# Name\nEntity type name, e.g., "Molecule".',
      '# Matching expression\nDefines objects in this type, e.g., "semtype=molecule".'
    ]);

    await this.action(
      'Explore entity type dialog',
      fromEvent(doneBtn, 'click'),
      undefined,
      'Click "Next" to proceed.'
    );

    // --- Step 5: Fill in entity type fields ---
    await this.dlgInputAction(typeDialog, 'Set "Name" to "molecule-tutorial"', 'Name', 'molecule-tutorial');
    await this.dlgInputAction(typeDialog, 'Set "Matching expression" to "semtype=Molecule"', 'Matching expression', 'semtype=Molecule');

    this.describe('Click OK to save the entity type.');
    await this.action('Save entity type', typeDialog.onClose, $(typeDialog.root).find('button.ui-btn-ok')[0]);

    // --- Step 6: Open Schemas Node ---
    this.title('Explore schemas');
    const schemasNode = stickyMetaNode ? findTreeNode(stickyMetaNode, 'Schemas', true) : undefined;
    this.describe('Schemas define the metadata fields and are linked to entity types. ' +
      'They are located under Browse → Platform → Sticky Meta → Schemas.');
    await this.action(
      'Open schemas node',
      waitForElementClick(schemasNode!.captionLabel),
      schemasNode!.captionLabel
    );
    await new Promise(resolve => setTimeout(resolve, 500));

    // --- Step 7: Open new schema dialog ---
    this.describe('Click "New Schema..." to create a schema.');
    const newSchemaButton = Array.from(document.querySelectorAll('button'))
      .find(b => b.textContent?.trim() === 'New Schema...');
    const schemaDialog = await this.openDialog(
      'Create a new schema',
      'Create a new schema',
      newSchemaButton,
      'Click "New Schema..." to start'
    );

    // --- Step 8: Explore schema dialog ---
    this.describe('Fill in the schema fields.');
    const schemaDialogRoot = schemaDialog.root;
    const schemaNameInput = Array.from(schemaDialogRoot.querySelectorAll('.ui-input-root'))
      .find(el => el.querySelector('label.ui-label span')?.textContent?.trim() === 'Name')
      ?.querySelector('input') as HTMLInputElement;
    const assocWithLabel = Array.from(schemaDialogRoot.querySelectorAll('.ui-input-root'))
      .find(el => el.querySelector('label.ui-label span')?.textContent?.trim() === 'Associated with:')
      ?.querySelector('div.ui-input-editor') as HTMLDivElement;
    const propertyTable = schemaDialogRoot.querySelector('table.d4-item-table') as HTMLTableElement;
    const propertyRow = propertyTable?.querySelectorAll('tbody > tr')[1] as HTMLTableRowElement; 

    doneBtn = describeElements([schemaNameInput, assocWithLabel, propertyRow], [
      '# Name\nName of the schema.',
      '# Associated with\nThe entity type this schema applies to.',
      '# Properties\nMetadata fields with a name and type (string, int, bool, double, datetime).'
    ]);

    await this.action('Explore schema dialog', fromEvent(doneBtn, 'click'));

    // --- Step 9: Fill in schema details ---
    await this.textInpAction(schemaDialogRoot, 'Set "Name" to "schema for tutorial"', 'Name', 'schema for tutorial');

    const assocWithRoot = $(schemaDialogRoot)
      .find('label.ui-label.ui-input-label span')
      .filter((_, el) => el.textContent?.trim() === 'Associated with:')
      .closest('.ui-input-root')[0]!;

    const selectEntitiesLabel = assocWithRoot.querySelector('.d4-link-action') as HTMLElement;
    await this.action(
      'Select associated entity',
      new Promise<void>((resolve) => selectEntitiesLabel.addEventListener('click', () => resolve(), { once: true })),
      selectEntitiesLabel
    );

    const selectorDlg = await new Promise<DG.Dialog>((resolve) => {
      const iv = setInterval(() => {
        const dlg = DG.Dialog.getOpenDialogs().find(d => d.title.includes('Select types for'));
        if (dlg) { clearInterval(iv); resolve(dlg); }
      }, 50);
    });

    const moleculeItem = Array.from(selectorDlg.root.querySelectorAll('.property-grid-item-name-text'))
      .find(x => x.textContent?.trim() === 'molecule-tutorial');
    const checkbox = moleculeItem?.closest('.property-grid-item')?.querySelector('input[type="checkbox"]') as HTMLInputElement;

    await this.action(
      'Select molecule-tutorial',
      new Promise<void>((resolve) => checkbox.addEventListener('click', () => resolve(), { once: true })),
      checkbox
    );

    const selectorOkBtn = $(selectorDlg.root).find('button.ui-btn-ok')[0];
    await this.action('Confirm entity selection', selectorDlg.onClose, selectorOkBtn);

    await this.textInpAction(propertyRow, 'Set property "Name" to "project name"', 'Name', 'project name');
    await this.choiceInputAction(propertyRow, 'Set property "Type" to "string"', 'Property Type', 'string');

    this.describe('Click OK to save schema.');
    await this.action('Save schema', schemaDialog.onClose, $(schemaDialog.root).find('button.ui-btn-ok')[0]);

    // --- Step 10: Open dataset ---
    this.title('Annotate dataset');
    this.t = await grok.data.loadTable(`${_package.webRoot}files/smiles.csv`);
    const tv = grok.shell.addTableView(this.t);
    this.describe('Blue circles indicate metadata availability. Click a cell to annotate.');
    grok.shell.windows.showContextPanel = true;

    const grid = tv.grid;
    const rowIndex = 0;
    const colName = 'smiles';
    await new Promise<void>((resolve) => {
      const sub = grid.onCellClick.subscribe((args) => {
        if (args.cell.rowIndex === rowIndex && args.cell.column.name === colName) { sub.unsubscribe(); resolve(); }
      });
    });

    await new Promise(resolve => setTimeout(resolve, 500));

    // --- Step 11: Fill in sticky meta property ---
    this.describe('Enter values for schema properties in the context panel.');
    const stickyPane = document.querySelector('.d4-accordion-pane[d4-title="Sticky meta"]') as HTMLElement;
    const header = stickyPane.querySelector('.d4-accordion-pane-header') as HTMLElement;
    if (header && !stickyPane.children[0].classList.contains('expanded')) header.click();

    const schemaSection = Array.from(stickyPane.querySelectorAll('.grok-sticky-meta-schema-header'))
      .find(el => el.textContent?.trim() === 'Schema for tutorial') as HTMLElement;

    const projectPropertyInput = schemaSection.parentElement?.querySelector('input[name="input-Project-name"]') as HTMLInputElement;
    await this.action(
      'Enter value for project name',
      new Promise<void>((resolve) => {
        const listener = () => {
          if (projectPropertyInput.value.trim() !== '') { 
            projectPropertyInput.removeEventListener('input', listener); resolve();
          }
        };
        projectPropertyInput.addEventListener('input', listener);
      }),
      projectPropertyInput
    );

    // --- Step 12: Save sticky meta ---
    const accordionPane = schemaSection.closest<HTMLDivElement>('.d4-accordion-pane');
    const saveBtn = Array.from(accordionPane!.querySelectorAll<HTMLButtonElement>('button[name^="button-Save"]'))
      .find(btn => !btn.classList.contains('disabled'));

    await this.action(
      'Save sticky meta changes',
      new Promise<void>((resolve) => saveBtn?.addEventListener('click', () => resolve(), { once: true })),
      saveBtn
    );

    // --- Step 13: Hover to see tooltip ---
    this.describe('Hover a cell to verify metadata tooltip.');
    await new Promise<void>((resolve) => {
      const sub = grid.onCellTooltip((cell) => {
        if (cell.cell.column.name === 'smiles' && cell.gridRow === 0) { sub.unsubscribe(); resolve(); }
      });
    });

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