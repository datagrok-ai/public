/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}



//name: testDep
export async function testDeploymnet(): Promise<null> {
  console.log("This package deployed successfully.");
  return null

}

// //output: viewer result
// export async function useViwer(): Promise<DG.DataFrame> {

//   const fastaDF = (await (grok.functions.eval(`OpenServerFile("System:AppData/Bio/samples/FASTA.fasta")`)))[0];
//   let view = grok.shell.addTableView(fastaDF);
//   const grid = view.grid;
//   for (var col in grid.columns){
//     console.log("Column:", col);
//   }

//   return fastaDF;
// }

//name: useViewer
//output: viewer result
export async function useViewer(): Promise<any> {
  // Open the FASTA file
  const fastaDF = (await (grok.functions.eval(`OpenServerFile("System:AppData/Bio/samples/FASTA.fasta")`)))[0];
  let view = grok.shell.addTableView(fastaDF);
  const grid = view.grid;
  
  // Find the macromolecule (sequence) column
  let sequenceColumn = null;
  let maxSequenceLength = 0;
  
  // Correctly iterate through columns
  for (let i = 0; i < grid.columns.length; i++) {
    const column = grid.columns.byIndex(i);
    if (column && column.column) {
      console.log("Column:", column.column.name, "SemType:", column.column.semType, "Tags:", column.column.tags);
      
      // Check if this is a macromolecule column (sequence column)
      if (column.column.semType === 'Macromolecule') {
        sequenceColumn = column;
        
        // Find max sequence length - needed for MSA header
        const values = fastaDF.getCol(column.column.name);
        for (let j = 0; j < values.length; j++) {
          const sequence = values[j];
          if (sequence && sequence.length > maxSequenceLength) {
            maxSequenceLength = sequence.length;
          }
        }
      }
    }
  }
  
  if (sequenceColumn) {
    // @ts-ignore
    console.log("Found sequence column:", sequenceColumn.column.name);
    console.log("Max sequence length:", maxSequenceLength);
    
    // Create MSA header for this column
    // addMSAHeaderToSequenceColumn(grid, sequenceColumn, maxSequenceLength);
  } else {
    console.error("No macromolecule/sequence column found in this FASTA file");
  }
  
  return view;
}