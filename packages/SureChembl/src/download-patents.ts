import * as grok from 'datagrok-api/grok';
import {Abstract, ClaimResponse, Description, PatentDocument, SureChemblBatchResponse} from './surechembl-api';
import * as DG from 'datagrok-api/dg';

const BATCH_SIZE = 1;
const SURE_CHEMBL_API_URL = 'https://www.surechembl.org/api';
export const S3_BASE_PATH = 'System:AppData/SureChembl/Patents/';

export interface PatentContentFile {
  document: PatentDocument;
  content: Uint8Array | string;
  isPdf: boolean;
}

export interface DownloadPatentStatus {
    id: string;
    status: string;
}

export async function downloadPatentDocuments(docIds: string[]): Promise<DownloadPatentStatus[]> {
  const progressBar = DG.TaskBarProgressIndicator.create(`Starting patents download...`);
  console.log(`~~~~~~~~~~~started dowloading patents`);
  const t = Date.now();
  const downloadStatus: DownloadPatentStatus[] = [];
  // Process documents in batches of 100
  for (let i = 0; i < docIds.length; i += BATCH_SIZE) {
    const batch = docIds.slice(i, i + BATCH_SIZE);
    const batchIds = batch.join(',');
    try {
      // Fetch patent documents
      const time1 = Date.now();
      const response = await grok.dapi.fetchProxy(`${SURE_CHEMBL_API_URL}/document/batch?doc_ids=${batchIds}`, {
        method: 'POST',
      });
      console.log(`${i} ******** received data for ${batch.length} patents for ${Date.now() - time1} ms`);
      const data: SureChemblBatchResponse = await response.json();

      if (data.status !== 'OK') {
        downloadStatus.push({id: batchIds, status: `Surechembl returned error: ${data.error_message}`});
        continue;
      }

      // Process each document in the batch
      for (const document of data.data) {
        const docId = document.doc_id;
        const s3Path = `${S3_BASE_PATH}${docId}`;

        // Try to fetch PDF if available
        if (document.contents.patentDocument.hasPDF) {
          try {
            const time2 = Date.now();
            const pdfResponse = await grok.dapi.fetchProxy(
              `${SURE_CHEMBL_API_URL}${document.contents.patentDocument.pdfUrl}`,
            );
            console.log(`${i} ******** received pdf for ${document.doc_id} for ${Date.now() - time2} ms`);
            const pdfData = await pdfResponse.arrayBuffer();
            const time3 = Date.now();
            await grok.dapi.files.write(`${s3Path}.pdf`, Array.from(new Uint8Array(pdfData)));
            downloadStatus
              .push({id: batchIds, status: `${s3Path}.pdf downloaded successfully in ${Date.now() - t} ms`});
            console.log(`${i} ******** wrote pdf for ${document.doc_id} for ${Date.now() - time3} ms`);
          } catch (error: any) {
            // If PDF fetch fails, fall back to text content
            const textContent = getFormattedPatentText(document);
            const time4 = Date.now();
            await grok.dapi.files.writeAsText(`${s3Path}.txt`, textContent);
            downloadStatus
              .push({id: batchIds, status: `${s3Path}.txt downloaded successfully in ${Date.now() - t} ms`});
            console.log(`${i} ******** PDF failed!!!, error: ${error.message ?? error}
                wrote text for ${document.doc_id} for ${Date.now() - time4} ms`);
          }
        } else {
          // If PDF is not available, use text content
          const textContent = getFormattedPatentText(document);
          const time5 = Date.now();
          await grok.dapi.files.writeAsText(`${s3Path}.txt`, textContent);
          downloadStatus
            .push({id: batchIds, status: `${s3Path}.txt downloaded successfully in ${Date.now() - t} ms`});
          console.log(`${i} ******** PDF failed!!! wrote text for ${document.doc_id} for ${Date.now() - time5} ms`);
        }
      }
    } catch (e: any) {
      downloadStatus.push({id: batchIds, status: `Error: ${e.message ?? e}`});
    } finally {
      progressBar?.update((i + 1)/docIds.length*100, `Downloaded ${i + 1} patents of ${docIds.length}`);
    }
  }
  console.log(`~~~~~~~~~~~downloaded ${docIds.length} patents for ${Date.now() - t}`);
  progressBar.close();
  return downloadStatus;
}

function getFormattedPatentText(
  document: PatentDocument,
  cleanHTML: boolean = true,
): string {
  const title = document.contents.patentDocument.bibliographicData.technicalData
    .inventionTitles?.[0]?.title || 'Untitled Patent';

  const doc = document.contents.patentDocument;

  const getContatenatedContent = (content: Abstract[] | Description[] | ClaimResponse[]) => {
    let contentStr = '';
    for (const item of content) {
      if (item.section?.content)
        contentStr += `${item.section?.content}'\n\n'`;
    }
    return contentStr;
  };

  const sections = [
    {name: 'TITLE', content: title},
    {name: 'ABSTRACT', content: getContatenatedContent(doc.abstracts ?? [])},
    {name: 'DESCRIPTION', content: getContatenatedContent(doc.descriptions ?? [])},
    {name: 'CLAIMS', content: getContatenatedContent(doc.claimResponses ?? [])},
  ];

  return sections
    .filter((s) => s.content)
    .map((s) => {
      let content = s.content;
      if (cleanHTML)
        content = content.replace(/<\/?[^>]+(>|$)/g, '');
      return `=== ${s.name} ===\n${content}`;
    })
    .join('\n\n');
}
