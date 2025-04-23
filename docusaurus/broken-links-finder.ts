import * as fs from 'fs';
import * as path from 'path';
import * as https from 'https';
import { argv } from 'process';

//Example of run:
//tsc broken-links-finder.ts
//node broken-links-finder.js 'C:\\Users\\sssav\\source\\repos\\public\\help'
const rootDir = argv[2];
const httpLink = /https:\/\/[^\s"'<>\)\]]+/g;
const isLinkValidMap = new Map<string, boolean>();
const linksInFile = new Map<string, string[]>();
const TIMEOUT_MS = 10000;

async function readFilesRecursively(dir: string): Promise<void> {
  try {
    const entries = await fs.promises.readdir(dir, { withFileTypes: true });
    for (const entry of entries) {
      const fullPath = path.join(dir, entry.name);

      if (entry.isDirectory()) {
        await readFilesRecursively(fullPath);
      } else if (entry.isFile() && (fullPath.endsWith('.md') || fullPath.endsWith('.mdx'))) {

        try {
          const data = await fs.promises.readFile(fullPath, 'utf8');
          const rawLinks = data.match(httpLink) || [];

          const cleanLinks = Array.from(
            new Set(rawLinks.map(link => link.split('#')[0].split('?')[0]))
          );

          for (const link of cleanLinks) {
            let status = await checkLinkExistence(link);
            isLinkValidMap.set(link, status);
            if (!status) {
              const current = linksInFile.get(fullPath) || [];
              current.push(link);
              linksInFile.set(fullPath, current);
            }
          }
        } catch (fileErr) {
          console.error(`Failed to read file ${fullPath}:`, fileErr);
        }
      }
    }
  } catch (dirErr) {
    console.error(`Failed to read directory ${dir}:`, dirErr);
  }
}

async function checkLinkExistence(url: string): Promise<boolean> {
  if (isLinkValidMap.has(url))
    return isLinkValidMap.get(url)!;

  return new Promise(async (resolve) => {
    try {
      const req = https.get(url, (res) => {
        const exists = res.statusCode !== undefined && res.statusCode < 400;
        resolve(exists);
      });

      req.setTimeout(TIMEOUT_MS, () => {
        req.destroy();
      });

      req.on('error', (e) => { console.log(`error on urlL: ${url} ${e}`); resolve(false) });
    } catch (err) {
      console.log(`catch ${url}`);
      resolve(false);
    }
  });
}

(async () => {
  await readFilesRecursively(rootDir);
  const obj = Object.fromEntries(linksInFile.entries());
  fs.writeFileSync('broken-links.json', JSON.stringify(obj, null, 2), 'utf8');
})();