import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { _package } from '../package';

export function showHelpPanel(info: string) {
  grok.shell.windows.help.visible = true;
  grok.shell.windows.help.showHelp(ui.markdown(info));
    
  grok.shell.windows.context.visible = true;  
  grok.shell.windows.showContextPanel = true;
  grok.shell.windows.showProperties = true;
  grok.shell.windows.help.visible = true;
}

export async function openDataset(name: string): Promise<DG.TableView> {
  const df = await _package.files.readAsText(name);
  const table = DG.DataFrame.fromCsv(df);
  return grok.shell.addTableView(table);
}

export const ESMFOLD_HELP = `# Folding Demo

The Folding Demo showcases EsmFold, a model that predicts the 3D structure of proteins from their amino acid sequences.

You can use EsmFold in two ways:

1. **For the entire column**:
   - Navigate to the top menu and select Bio | BioNemo | EsmFold.
   - This will apply EsmFold to all sequences in the selected column.

2. **For a single cell**:
   - Click on a specific sequence in the dataset.
   - The ESMFold panel will appear in the property panel.`

export const DIFFDOCK_HELP = `# DiffDock Demo

The DiffDock Demo showcases DiffDock, a model that predicts the 3D structure of how a molecule interacts with a protein.

You can use DiffDock to explore molecular interactions in two ways:

1. **For viewing poses**:
   - Click on a pose in the dataset.
   - The Mol* viewer will open, displaying the selected pose.

2. **For viewing other poses and confidence**:
   - In the Mol* viewer, a combo popup will list other generated poses and their confidence levels.
   - To visualize other generated poses, click on the pose of your interest.`