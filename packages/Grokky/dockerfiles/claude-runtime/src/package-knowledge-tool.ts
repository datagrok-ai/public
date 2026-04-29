import * as fs from 'node:fs/promises';
import * as path from 'node:path';
import * as YAML from 'yaml';
import {z} from 'zod/v4';
import {createSdkMcpServer, tool} from '@anthropic-ai/claude-agent-sdk';
import {WORKSPACE} from './sync-utils';

export interface PackageKnowledge {
  packageName: string;
  description: string;
  keywords: string[];
  overview?: string;
  apiRef?: string;
  docsRef?: string;
}

export interface LoadedPackage extends PackageKnowledge {
  apiRefPath?: string;
  docsRefPath?: string;
}

const KNOWLEDGE_FILE = 'package-knowledge.yaml';

let cache: Map<string, LoadedPackage> | null = null;

export function invalidatePackageKnowledgeCache(): void {
  cache = null;
}

async function readKnowledge(packagesDir: string, dir: string): Promise<LoadedPackage | null> {
  const knowledgePath = path.join(packagesDir, dir, 'agents', KNOWLEDGE_FILE);
  try {
    const raw = await fs.readFile(knowledgePath, 'utf-8');
    const parsed = YAML.parse(raw) as PackageKnowledge;
    if (!parsed?.packageName || !parsed.description) {
      console.warn(`package-knowledge: ${dir}/${KNOWLEDGE_FILE} missing required fields, skipping`);
      return null;
    }
    const base = path.join(packagesDir, dir);
    return {
      ...parsed,
      apiRefPath: parsed.apiRef ? path.resolve(base, parsed.apiRef) : undefined,
      docsRefPath: parsed.docsRef ? path.resolve(base, parsed.docsRef) : undefined,
    };
  } catch { /* no knowledge file for this package */ }
  return null;
}

export async function loadPackageKnowledge(): Promise<Map<string, LoadedPackage>> {
  if (cache)
    return cache;
  const map = new Map<string, LoadedPackage>();
  const packagesDir = path.join(WORKSPACE, 'packages');
  let pkgDirs: import('node:fs').Dirent[];
  try {
    pkgDirs = await fs.readdir(packagesDir, {withFileTypes: true});
  } catch (e: any) {
    console.warn(`package-knowledge: failed to scan ${packagesDir}: ${e.message}`);
    cache = map;
    return map;
  }
  const loaded = await Promise.all(
    pkgDirs.filter((e) => e.isDirectory()).map((e) => readKnowledge(packagesDir, e.name)),
  );
  for (const pkg of loaded) {
    if (pkg)
      map.set(pkg.packageName, pkg);
  }
  cache = map;
  return map;
}

const getPackageKnowledgeTool = tool(
  'get_package_knowledge',
  'Returns pointers to the authoritative reference files for a Datagrok package: ' +
  'description, keywords, and absolute workspace paths to the API reference ' +
  '(function signatures + types) and docs. Call this BEFORE exploring with Glob/Grep. ' +
  'Then use Read on the returned apiRef path to locate the function you need — do NOT ' +
  'try to discover package APIs by browsing the filesystem. The packageName must match ' +
  'an entry from the "## Available Packages" table in the system prompt (case-sensitive).',
  {packageName: z.string().describe('Exact package name from the index (e.g. "Chem", "Bio")')},
  async ({packageName}) => {
    const map = await loadPackageKnowledge();
    const pkg = map.get(packageName);
    if (!pkg) {
      const known = [...map.keys()].sort().join(', ');
      return {content: [{type: 'text' as const, text: `Unknown package "${packageName}". Known packages: ${known}`}]};
    }
    const text =
      `# ${pkg.packageName}\n\n${pkg.description}\n\n` +
      (pkg.overview ? `## Overview\n\n${pkg.overview}\n\n` : '') +
      `Keywords: ${pkg.keywords.join(', ')}\n\n` +
      (pkg.apiRefPath ? `apiRef: ${pkg.apiRefPath}\n` : '') +
      (pkg.docsRefPath ? `docsRef: ${pkg.docsRefPath}\n` : '') +
      (pkg.apiRefPath || pkg.docsRefPath ? `\nNext step: Read the apiRef path to find the function you need.\n` : '');
    return {content: [{type: 'text' as const, text}]};
  },
);

// Create a fresh SDK MCP server instance per call — McpServer's Protocol
// allows only one transport at a time, and the SDK calls connect() on every query().
export function createPackageKnowledgeServer() {
  return createSdkMcpServer({
    name: 'datagrok-knowledge',
    version: '1.0.0',
    tools: [getPackageKnowledgeTool],
  });
}
