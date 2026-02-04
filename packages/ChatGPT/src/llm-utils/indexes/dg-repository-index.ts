/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {LLMClient} from '../LLM-client';
import OpenAI from 'openai';

const VECTOR_STORE_NAME = 'datagrok-public-index';
const MAX_CONCURRENCY = 3;
const MAX_RETRIES = 5;
const BASE_RETRY_DELAY_MS = 1000;

function sleep(ms: number): Promise<void> {
  return new Promise((resolve) => setTimeout(resolve, ms));
}

function isRetryableError(err: any): boolean {
  const status = err?.status ?? err?.response?.status ?? err?.statusCode;
  if (status === 429)
    return true;
  if (status >= 500 && status < 600)
    return true;
  const code = err?.code ?? err?.cause?.code;
  return code === 'ETIMEDOUT' || code === 'ECONNRESET' || code === 'ENOTFOUND';
}

async function withRetry<T>(fn: () => Promise<T>, label: string): Promise<T> {
  let attempt = 0;
  while (true) {
    try {
      return await fn();
    } catch (err) {
      attempt++;
      if (attempt >= MAX_RETRIES || !isRetryableError(err))
        throw err;
      const delay = BASE_RETRY_DELAY_MS * Math.pow(2, attempt - 1) + Math.floor(Math.random() * 250);
      console.warn(`Retrying ${label} (attempt ${attempt}/${MAX_RETRIES}) after ${delay}ms`, err);
      await sleep(delay);
    }
  }
}

async function getOrCreateVectorStoreId(client: OpenAI, name: string, pg: DG.TaskBarProgressIndicator): Promise<string> {
  console.log(`Checking for existing vector store with name: ${name}`);
  pg.update(0, `Checking for existing vector store`);
  let page: OpenAI.VectorStores.VectorStoresPage | null = null;
  while (true) {
    if (!page)
      page = await client.vectorStores.list();
    else
      page = await page.getNextPage();
    const found = page.data?.find((vs) => vs.name === name);
    if (found) {
      console.log(`Found existing vector store with name: ${name}, id: ${found.id}`);
      return found.id;
    }
    if (!page.has_more)
      break;
  }
  console.log(`No existing vector store found with name: ${name}. Creating a new one.`);
  const created = await client.vectorStores
    .create({name, description: 'Datagrok public repository vector store with all indexable files (code + docs).'});
  console.log(`Created new vector store with name: ${name}, id: ${created.id}`);
  DG.Utils.download('vector-store-id.txt', created.id);
  return created.id;
}

async function clearVectorStore(client: OpenAI, vectorStoreId: string, pg: DG.TaskBarProgressIndicator): Promise<void> {
  let page: OpenAI.VectorStores.Files.VectorStoreFilesPage | undefined = undefined;
  console.log(`Clearing vector store with id: ${vectorStoreId}`);
  pg.update(0, `Clearing vector store`);
  const processedSet = new Set<string>();
  while (true) {
    if (!page)
      page = await client.vectorStores.files.list(vectorStoreId, {limit: 100});
    else
      page = await page.getNextPage();
    const files = page.data ?? [];
    const count = files.length;
    let loopCount = 0;
    const curProcessedSize = processedSet.size;
    for (const f of files) {
      if (processedSet.has(f.id))
        continue;
      processedSet.add(f.id);
      pg.update(loopCount / count * 100, `Deleting file ${loopCount + 1} of ${count}`);
      try {
        await client.files.delete(f.id); // delete file and also remove it from any vector stores.
      } catch (e) {
        console.error(`Error deleting file : ${f.id}`, e);
      }
      loopCount++;
    }
    if (processedSet.size === curProcessedSize) {
      console.log('No more files to delete.');
      // openai API does some weird thing with caching.... and has more flag is not reliable
      break;
    }
    if (!page.has_more)
      break;
  }
}

export async function listIndexableFiles(owner: string, repo: string): Promise<string[]> {
  const repoInfo = await grok.dapi.fetchProxy(
    `https://api.github.com/repos/${owner}/${repo}`
  ).then((r) => r.json());
  const sha = repoInfo.default_branch;

  const tree = await grok.dapi.fetchProxy(
    `https://api.github.com/repos/${owner}/${repo}/git/trees/${sha}?recursive=1`
  ).then((r) => r.json());

  const allowedExtensions = [
    '.ts',
    '.tsx',
    '.md',
    '.mdx',
    '.js',
    '.py',
    '.sql',
    '.r',
    '.cpp',
    '.h',
  ];
  const isAllowedExtension = (path: string) =>
    allowedExtensions.some((ext) => path.endsWith(ext));
  const nonAllowedDirs = [
    '/node_modules/',
    '/dist/',
    '/build/',
    '.git/',
    '/.github/',
    'grok_connect',
    'docusaurus',
  ];
  const isInNonAllowedDir = (path: string) =>
    nonAllowedDirs.some((dir) => path.includes(dir));
  return tree.tree
    .filter(
      (item: any) =>
        item.type === 'blob' &&
        isAllowedExtension(item.path) &&
        !isInNonAllowedDir(item.path)
    )
    .map((item: any) => item.path);
}

export async function uploadFilesToVectorStroreOneByOne() {
  const pg = DG.TaskBarProgressIndicator.create(
    'Uploading files to vector store...'
  );
  const links = await listIndexableFiles('datagrok-ai', 'public').then(
    (res) => {
      // save to a file
      /**  */
      const downloadLinks = res.map(
        (p) =>
          `https://raw.githubusercontent.com/datagrok-ai/public/master/${p}`
      );
      DG.Utils.download('indexed_files.txt', downloadLinks.join('\n'));
      //   console.log('indexed_files.txt', downloadLinks.join('\n'));
      return downloadLinks;
    }
  );
  // go one by one and upload
  const client = LLMClient.getInstance().openai;

  // check if the vector store with name 'datagrok-public-index' already exists, if not, create it
  const vectorStoreId = await getOrCreateVectorStoreId(client, VECTOR_STORE_NAME, pg);
  await clearVectorStore(client, vectorStoreId, pg);

  // const file = await client.vectorStores.files.uploadAndPoll(
  //     process.env.VECTOR_STORE_ID,
  //     fs.createReadStream('C:/Users/rizhi/Desktop/GROK/public/help/learn/custom-machine-learning-models.md')
  // );

  const errors: any[] = [];
  let count = 0;
  const processedIds: string[] = [];
  let nextIndex = 0;

  const processLink = async (link: string) => {
    const response = await withRetry(() => fetch(link), 'download');
    if (!response || !response.ok) {
      throw new Error(
        `Failed to download file: ${response ? response.statusText : 'No response'}`
      );
    }

    const buffer = await response.arrayBuffer();
    const rawBase = 'https://raw.githubusercontent.com/datagrok-ai/public/master/';
    const repoPath = link.substring(rawBase.length);
    const normalizedPath = repoPath
      .replace('.mdx', '.md')
      .replace('.tsx', '.ts')
      .replace('.sql', '.txt')
      .replace('.r', '.txt')
      .replace('.h', '.c')
      .replace('.ts', '.js')// )))) apparently they made ts unsupported in vector store... so we upload as js
      .replace('.cpp', '.c')
      .replace('.c', '.txt');

    const fileName = normalizedPath;
    const originalExtension = repoPath.includes('.') ?
      repoPath.split('.').pop()! :
      'txt';
    const lastSlash = normalizedPath.lastIndexOf('/');
    const dir = lastSlash >= 0 ? normalizedPath.slice(0, lastSlash) : '';

    const file = new File([buffer], fileName);

    const fileObj = await withRetry(() => client.files.create({
      file,
      purpose: 'assistants',
    }), 'files.create');

    const fileRes = await withRetry(() => client.vectorStores.files.create(vectorStoreId, {
      file_id: fileObj.id, attributes: {
        repoPath: repoPath,
        originalExtension: originalExtension,
        directory: dir,
        firstParentFolder: dir.split('/')[0] ?? '', // good for JS-api specific searches
        secondParentFolder: dir.split('/')[1] ?? '', // good for libraries / packages
        originalFileName: repoPath.includes('/') ?
          repoPath.slice(repoPath.lastIndexOf('/') + 1) :
          repoPath,
      }
    }), 'vectorStores.files.create');
    processedIds.push(fileRes.id);
  };

  const workers = Array.from({length: Math.min(MAX_CONCURRENCY, links.length)}, async () => {
    while (true) {
      const i = nextIndex++;
      if (i >= links.length)
        break;
      const link = links[i];
      try {
        await processLink(link);
      } catch (err) {
        console.error('Error uploading file from link:', link, err);
        errors.push({link, error: err ? err : 'unknown error'});
      }
      count++;
      console.log(`progressed so far: ${count}/${links.length}`);
      pg.update((count / links.length) * 100, `Uploaded ${count} of ${links.length} files`);
    }
  });

  await Promise.all(workers);

  pg.close();


  console.log(`\n✓ Uploaded ${processedIds.length} of ${links.length} files`);

  // save errors to a file
  DG.Utils.download('upload_errors.json', JSON.stringify(errors, null, 2));
  // save processed ids to a file
  DG.Utils.download('processed_file_ids.txt', processedIds.join('\n'));
  // download the vector store id
  grok.shell.info(`File upload to vector store completed with ${errors.length} / ${processedIds.length} errors.`);
  grok.shell.info('Make sure to update the vector store id in the CHATGPT plugin settings and refresh the page.');
}
