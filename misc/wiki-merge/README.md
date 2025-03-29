# Documentation Copier

A TypeScript command-line utility to copy GitHub project documentation (markdown files and images) to a new location, while only including files that are referenced from top-level markdown files.

## Features

- Copies only documentation files (.md, .gif, .png) that are referenced in top-level markdown files
- Follows references in markdown files (links, images)
- Renames README.md files to match their containing folder name
- Preserves directory structure
- Supports processing multiple projects with a single command
- References TypeScript files via GitHub URLs instead of copying them
- Respects source directory boundaries (doesn't process files outside specified source directory)
- Organizes output by category based on package.json
- Simple command-line interface

## Installation

1. Ensure you have Node.js and TypeScript installed
2. Compile the TypeScript file:
   ```bash
   tsc wiki-merge.ts --target ES2018 --module commonjs --esModuleInterop
   ```

## Usage

### Single Project Mode

```bash
node wiki-merge.js <source-folder> <destination-folder>
```

Example:
```bash
node wiki-merge.js ./my-project/docs ./documentation-site/content
```

### Multiple Projects Mode

Process all projects under a root folder:

```bash
node wiki-merge.js --all <root-source-folder> <root-destination-folder>
```

Example:
```bash
node wiki-merge.js --all ./github-projects ./documentation-site/projects
```

In this mode, the script will:
1. Find all directories immediately under the specified root source folder
2. Process each directory as a separate project (skipping ones without README.md)
3. Maintain the same directory structure in the destination folder
4. Organize files by category (read from each project's package.json)
5. Provide a summary of successfully processed projects

## How It Works

The script:
1. Identifies all top-level markdown files in the source directory
2. Parses these files to find all references to:
   - Other markdown files via links: `[text](link.md)`
   - Images via markdown syntax: `![alt](image.png)`
   - Images via HTML tags: `<img src="image.gif" />`
3. Also finds references in any referenced markdown files (one level deep)
4. Copies only the top-level markdown files and their referenced files
5. Renames README.md files to match their containing folder name
6. References to TypeScript files (.ts) are converted to GitHub URLs using the format: 
   `https://github.com/datagrok-ai/public/tree/master/packages/{packageName}/{path-to-file}`
7. Skips files that are outside the specified source directory
8. Organizes output by category, reading the `category` property from each project's package.json
   - If the category property is missing, files are placed in an "Uncategorized" category

## Supported Reference Types

The script recognizes the following types of references in markdown files:

- Markdown links: `[Link text](path/to/file.md)`
- Markdown images: `![Alt text](path/to/image.png)`
- HTML img tags: `<img src="path/to/image.gif" />`
- TypeScript references: `[Link text](path/to/file.ts)` (converted to GitHub URLs)

## Category Organization

Projects are organized in the output directory based on their category:

```
output/
├── Chemistry/
│   ├── Chem.md
│   ├── docs/
│   │   └── api.md
│   └── images/
│       └── diagram.png
├── Geo/
│   ├── GIS.md
│   └── docs/
│       └── api.md
└── Uncategorized/
    └── MiscUtils.md
```

The category is determined by reading the `category` property in each project's package.json file:

```json
{
  "name": "chem",
  "version": "1.0.0",
  "description": "Chemistry package for Datagrok",
  "category": "Chemistry",
  "dependencies": {
    "datagrok-api": "^1.0.0"
  }
}
```

If the `category` property is missing or the package.json file doesn't exist, the project is placed in the "Uncategorized" category.

## Notes

- Only files referenced in top-level markdown files (and their direct references) will be copied
- External URLs (http://, https://, ftp://) are ignored
- Any non-documentation files (.txt, .js, etc.) will be ignored even if referenced
- TypeScript files (.ts) are not copied but instead referenced via GitHub URLs
- References to files outside the source directory are kept as-is in the markdown but will not be copied
- External directory references (e.g., with `../` or absolute paths) that point outside the source directory are skipped
- Existing files in the destination directory with the same name will be overwritten
- Warnings will be displayed for referenced files that cannot be found
- In multiple projects mode, directories without a README.md file will be skipped
- Circular references between markdown files are automatically detected and handled to prevent infinite loops 