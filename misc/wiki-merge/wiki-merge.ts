#!/usr/bin/env node

import * as fs from 'fs';
import * as fsPromises from 'fs/promises';
import * as path from 'path';

// GitHub base URL for references
const GITHUB_BASE_URL = 'https://github.com/datagrok-ai/public/tree/master/packages';

// Default category for packages without a category property
const DEFAULT_CATEGORY = 'Uncategorized';

// Add at the top with other constants
const packageCategories: Map<string, string> = new Map();

// Function to check if a file is a markdown file
function isMarkdownFile(fileName: string): boolean {
  return path.extname(fileName).toLowerCase() === '.md';
}

// Function to check if a file is an image file
function isImageFile(fileName: string): boolean {
  const extension = path.extname(fileName).toLowerCase();
  return extension === '.gif' || extension === '.png';
}

// Function to check if a file is a TypeScript or JavaScript file.
// This is used to identify files that should be referenced via GitHub URLs.
function isGitHubFile(fileName: string): boolean {
  const extension = path.extname(fileName).toLowerCase();
  return extension === '.ts' || extension === '.js';
}

// Function to check if a path is an external URL
function isExternalUrl(filePath: string): boolean {
  return filePath.startsWith('http://') || 
         filePath.startsWith('https://') || 
         filePath.startsWith('ftp://');
}

// Function to copy a file from source to destination
async function copyFile(source: string, destination: string): Promise<void> {
  await fsPromises.mkdir(path.dirname(destination), { recursive: true });
  await fsPromises.copyFile(source, destination);
  console.log(`Copied: ${source} → ${destination}`);
}

// Get the category from package.json
async function getPackageCategory(sourcePath: string): Promise<string> {
  try {
    const packageJsonPath = path.join(sourcePath, 'package.json');
    
    if (!fs.existsSync(packageJsonPath)) {
      return DEFAULT_CATEGORY;
    }
    
    const packageJsonContent = await fsPromises.readFile(packageJsonPath, 'utf8');
    const packageJson = JSON.parse(packageJsonContent);
    
    const category = packageJson.category || DEFAULT_CATEGORY;
    return category.replace(/\s+/g, '-');
  } catch (error) {
    console.warn(`Warning: Could not read package.json or parse category: ${error}`);
    return DEFAULT_CATEGORY;
  }
}

// Convert a local file path to a GitHub URL
function convertToGitHubUrl(filePath: string, sourceRoot: string): string {
  const relativePath = path.relative(sourceRoot, filePath);
  const packageName = relativePath.split(path.sep)[0]; // First directory is the package name
  
  // For TypeScript files under a package, create GitHub URL
  return `${GITHUB_BASE_URL}/${packageName}/${relativePath.split(path.sep).slice(1).join('/')}`;
}

// Helper function to extract references using a regex pattern
function extractReferencesWithPattern(
  content: string, 
  basePath: string, 
  regex: RegExp, 
  groupIndex: number,
  references: Set<string>,
  sourceDirRoot: string
): void {
  let match;
  while ((match = regex.exec(content)) !== null) {
    const link = match[groupIndex].split('#')[0]; // Remove anchor part if present
    if (link && !isExternalUrl(link)) {
      const absolutePath = path.resolve(basePath, link);
      
      references.add(absolutePath);
    }
  }
}

// Check if a path is within a specified directory
function isWithinDirectory(filePath: string, directory: string): boolean {
  // Normalize paths to handle different path formats
  const normalizedFilePath = path.normalize(filePath);
  const normalizedDirectory = path.normalize(directory);
  
  // Check if the file path starts with the directory path
  return normalizedFilePath.startsWith(normalizedDirectory);
}

// Function to extract all references from markdown content
function extractReferencesFromMarkdown(content: string, basePath: string, sourceDirRoot: string): Set<string> {
  const references = new Set<string>();
  
  // Find markdown links: [text](link)
  extractReferencesWithPattern(content, basePath, /\[([^\]]+)\]\(([^)]+)\)/g, 2, references, sourceDirRoot);
  
  // Find image references: ![alt](image-path)
  extractReferencesWithPattern(content, basePath, /!\[([^\]]*)\]\(([^)]+)\)/g, 2, references, sourceDirRoot);
  
  // Find HTML img tags: <img src="path" />
  extractReferencesWithPattern(content, basePath, /<img[^>]*src=["']([^"']+)["'][^>]*>/g, 1, references, sourceDirRoot);
  
  return references;
}

// Add a new function to transform README links
function transformReadmeLinks(link: string): string {
  // Check if it's a relative path to a README.md
  const readmeMatch = link.match(/\.\.\/([^\/]+)\/README\.md/);
  if (readmeMatch) {
    const packageName = readmeMatch[1];
    const category = packageCategories.get(packageName);
    if (category) {
      return link.replace(`${packageName}/README.md`, `${category}/${packageName}.md`);
    }
  }
  return link;
}

// Modify replaceTypeScriptReferences to also handle README transformations
async function replaceTypeScriptReferences(
  content: string, 
  basePath: string, 
  sourceRoot: string,
  sourceDirRoot: string
): Promise<string> {
  let updatedContent = content;
  
  // Replace markdown links: [text](link)
  const markdownLinkRegex = /\[([^\]]+)\]\(([^)]+)\)/g;
  updatedContent = updatedContent.replace(markdownLinkRegex, (match, text, link) => {
    if (isExternalUrl(link)) return match;
    
    // Transform README links first
    const transformedLink = transformReadmeLinks(link);
    
    // Then handle TypeScript/JavaScript files if needed
    if (isGitHubFile(transformedLink)) {
      // Handle both relative and absolute paths
      let absolutePath: string;
      if (transformedLink.startsWith('../../')) {
        // For paths starting with ../../, resolve relative to the sourceRoot
        absolutePath = path.resolve(sourceRoot, transformedLink.slice(6)); // Remove ../../ from the start
      } else {
        // Handle relative paths
        absolutePath = path.resolve(basePath, transformedLink.split('#')[0]);
      }
      
      const githubUrl = convertToGitHubUrl(absolutePath, sourceRoot);
      return `[${text}](${githubUrl})`;
    }
    
    console.log(`[${text}](${transformedLink})`);
    return `[${text}](${transformedLink})`;
  });
  
  // Replace image references: ![alt](image.ts) - not common but just in case
  const imageRegex = /!\[([^\]]*)\]\(([^)]+\.ts[^)]*)\)/g;
  updatedContent = updatedContent.replace(imageRegex, (match, alt, link) => {
    if (isExternalUrl(link)) return match;
    
    const absolutePath = path.resolve(basePath, link.split('#')[0]);

    const githubUrl = convertToGitHubUrl(absolutePath, sourceRoot);
    return `![${alt}](${githubUrl})`;
  });
  
  return updatedContent;
}

// Extract and process references from a markdown file
async function processReferencesFromFile(
  filePath: string, 
  filesToCopy: Set<string>,
  tsFilesReferenced: Set<string>,
  sourceRoot: string,
  sourceDirRoot: string,
  processedFiles: Set<string> = new Set(),
  depth: number = 0
): Promise<void> {
  // Skip if already processed to prevent infinite loops
  if (processedFiles.has(filePath)) {
    return;
  }
  
  // Check if the file is within the source directory
  if (!isWithinDirectory(filePath, sourceDirRoot)) {
    console.warn(`Warning: Skipping file outside source directory: ${filePath}`);
    return;
  }
  
  // Mark as processed
  processedFiles.add(filePath);
  
  try {
    const content = await fsPromises.readFile(filePath, 'utf-8');
    const dirName = path.dirname(filePath);
    const references = extractReferencesFromMarkdown(content, dirName, sourceDirRoot);
    
    for (const reference of references) {
      if (fs.existsSync(reference)) {
        // If it's a TypeScript file, add to tsFilesReferenced instead of filesToCopy
        if (isGitHubFile(reference)) {
          tsFilesReferenced.add(reference);
        } else {
          filesToCopy.add(reference);
          
          // If reference is a markdown file, process its references too
          if (isMarkdownFile(reference)) {
            await processReferencesFromFile(
              reference, 
              filesToCopy, 
              tsFilesReferenced,
              sourceRoot,
              sourceDirRoot,
              processedFiles, 
              depth + 1
            );
          }
        }
      } else {
        console.warn(`Warning: Referenced file not found: ${reference}`);
      }
    }
  } catch (error) {
    console.warn(`Warning: Could not read file: ${filePath}`);
  }
}

// Process a markdown file and update its content to use GitHub URLs for TypeScript files
async function processMarkdownContent(
  sourceFilePath: string, 
  destFilePath: string,
  sourceRoot: string,
  sourceDirRoot: string,
  tsFilesReferenced: Set<string>
): Promise<void> {
  try {
    // Read the markdown content
    let content = await fsPromises.readFile(sourceFilePath, 'utf-8');
    const dirName = path.dirname(sourceFilePath);
    
    // Replace ../../help/ with ../../../ in the content before extracting references
    content = content.replace(/\.\.\/\.\.\/help\//g, "../../../");

    // Replace TypeScript file references with GitHub URLs
    const updatedContent = await replaceTypeScriptReferences(content, dirName, sourceRoot, sourceDirRoot);
    
    // Write the updated content
    await fsPromises.mkdir(path.dirname(destFilePath), { recursive: true });
    await fsPromises.writeFile(destFilePath, updatedContent, 'utf-8');
    
    console.log(`Processed: ${sourceFilePath} → ${destFilePath}`);
  } catch (error) {
    console.warn(`Warning: Could not process markdown file: ${sourceFilePath}`);
  }
}

// Main function to copy documentation
async function copyDocumentation(source: string, destination: string, category?: string): Promise<void> {
  try {
    // Normalize source path
    source = path.normalize(source);
    
    // Check if source directory exists
    if (!fs.existsSync(source)) {
      console.error(`Error: Source directory '${source}' does not exist.`);
      process.exit(1);
    }

    // Always get category from package.json, ignore passed category parameter
    category = await getPackageCategory(source);
    
    // Create destination directory with category
    const categoryDestination = path.join(destination, category);
    if (!fs.existsSync(categoryDestination)) {
      await fsPromises.mkdir(categoryDestination, { recursive: true });
    }
    
    // Find README.md in the top-level directory
    const readmePath = path.join(source, 'README.md');
    if (!fs.existsSync(readmePath)) {
      console.error(`Error: README.md not found in source directory: ${source}`);
      process.exit(1);
    }
    
    // Set to track files to be copied (non-TypeScript files)
    const filesToCopy = new Set<string>();
    
    // Set to track TypeScript files referenced (for GitHub URLs)
    const tsFilesReferenced = new Set<string>();
    
    // Get the top-level source directory (packages directory)
    const sourceRoot = path.resolve(source, '..');
    
    // Use source directory as the root for path restriction
    const sourceDirRoot = source;
    
    // Process references from README.md (with empty processedFiles set)
    await processReferencesFromFile(readmePath, filesToCopy, tsFilesReferenced, sourceRoot, sourceDirRoot, new Set());
    
    // Add README to the files that need markdown processing
    const markdownFiles = new Set([readmePath]);
    for (const filePath of filesToCopy) {
      if (isMarkdownFile(filePath)) {
        markdownFiles.add(filePath);
      }
    }
    
    console.log(`Found ${filesToCopy.size} files to copy and ${tsFilesReferenced.size} TypeScript files to reference via GitHub`);
    
    // Process markdown files to update TypeScript references
    for (const markdownFile of markdownFiles) {
      const relativePath = path.relative(source, markdownFile);
      let destPath = path.join(categoryDestination, relativePath);
      
      // Handle README.md special case
      if (markdownFile === readmePath) {
        const sourceFolderName = path.basename(source);
        destPath = path.join(categoryDestination, `${sourceFolderName}.md`);
      }
      
      await processMarkdownContent(markdownFile, destPath, sourceRoot, sourceDirRoot, tsFilesReferenced);
      
      // Remove from filesToCopy since we've handled it separately
      filesToCopy.delete(markdownFile);
    }
    
    // Copy remaining non-markdown and non-TypeScript files
    for (const fileToCopy of filesToCopy) {
      // Skip TypeScript files
      if (isGitHubFile(fileToCopy)) {
        continue;
      }
      
      // Compute relative path from source directory
      const relativePath = path.relative(source, fileToCopy);
      if (!relativePath.startsWith('..')) { // Ensure file is within source directory
        const destPath = path.join(categoryDestination, relativePath);
        await copyFile(fileToCopy, destPath);
      }
    }
    
    // Store the package name and category mapping
    const packageName = path.basename(source);
    packageCategories.set(packageName, category || await getPackageCategory(source));
    
    console.log(`Documentation successfully copied from '${source}' to '${categoryDestination}'.`);
    console.log(`Category: ${category}`);
    if (tsFilesReferenced.size > 0) {
      console.log(`${tsFilesReferenced.size} TypeScript files were referenced via GitHub URLs instead of being copied.`);
    }
  } catch (error) {
    console.error('Error copying documentation:', error);
    process.exit(1);
  }
}

// Process multiple projects under a root folder
async function processAllProjects(rootSourceDir: string, rootDestDir: string): Promise<void> {
  try {
    // Normalize root source path
    rootSourceDir = path.normalize(rootSourceDir);
    
    // Check if source directory exists
    if (!fs.existsSync(rootSourceDir)) {
      console.error(`Error: Root source directory '${rootSourceDir}' does not exist.`);
      process.exit(1);
    }

    // Delete the root destination directory if it exists
    if (fs.existsSync(rootDestDir)) {
      await fsPromises.rm(rootDestDir, { recursive: true, force: true });
      console.log(`Deleted existing directory: ${rootDestDir}`);
    }

    // Create fresh root destination directory
    await fsPromises.mkdir(rootDestDir, { recursive: true });
    
    // Get all subdirectories of the root source directory
    const items = await fsPromises.readdir(rootSourceDir, { withFileTypes: true });
    const projectDirs = items
      .filter(item => item.isDirectory())
      .map(item => item.name);
    
    if (projectDirs.length === 0) {
      console.warn(`Warning: No project directories found in '${rootSourceDir}'.`);
      return;
    }
    
    console.log(`Found ${projectDirs.length} project directories to process.`);
    
    // Process each project directory
    let successCount = 0;
    let errorCount = 0;
    const categoryCounts: Record<string, number> = {};
    
    for (const projectDir of projectDirs) {
      const projectSourcePath = path.join(rootSourceDir, projectDir);
      
      // Check if README.md exists in the project directory
      const readmePath = path.join(projectSourcePath, 'README.md');
      if (!fs.existsSync(readmePath)) {
        console.warn(`Warning: Skipping '${projectDir}' - No README.md found.`);
        errorCount++;
        continue;
      }
      
      try {
        console.log(`Processing project: ${projectDir}`);
        
        // Get the category from package.json
        const category = await getPackageCategory(projectSourcePath);
        
        // Update category counts
        categoryCounts[category] = (categoryCounts[category] || 0) + 1;
        
        await copyDocumentation(projectSourcePath, rootDestDir, category);
        successCount++;
      } catch (error) {
        console.error(`Error processing project '${projectDir}': ${error}`);
        errorCount++;
      }
    }
    
    console.log(`\nSummary:`);
    console.log(`  Projects processed successfully: ${successCount}`);
    console.log(`  Projects with errors: ${errorCount}`);
    console.log(`  Total projects: ${projectDirs.length}`);
    console.log(`\nCategory breakdown:`);
    Object.entries(categoryCounts).forEach(([category, count]) => {
      console.log(`  ${category}: ${count} projects`);
    });
  } catch (error) {
    console.error('Error processing projects:', error);
    process.exit(1);
  }
}

// Parse command-line arguments
function main(): void {
  const args = process.argv.slice(2);
  
  // Check for --all flag
  const allModeIndex = args.indexOf('--all');
  const allMode = allModeIndex !== -1;
  
  // Remove the --all flag if present
  if (allMode) {
    args.splice(allModeIndex, 1);
  }
  
  if (args.length !== 2) {
    console.error(`Usage: 
  - Single project mode: npx ts-node wiki-merge.ts <source-folder> <destination-folder>
  - All projects mode: npx ts-node wiki-merge.ts --all <root-source-folder> <root-destination-folder>`);
    process.exit(1);
  }
  
  const [source, destination] = args.map(arg => path.resolve(arg));
  
  if (allMode) {
    processAllProjects(source, destination);
  } else {
    copyDocumentation(source, destination);
  }
}

// Export functions for testing
export {
  replaceTypeScriptReferences,
  isGitHubFile,
  convertToGitHubUrl,
  isWithinDirectory,
  transformReadmeLinks
};

// Only run main if this is the entry point
if (require.main === module) {
  main();
} 