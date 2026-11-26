# Files

Playlist: **Files**. Feature name of every clip: `How-to · Files`. Features in
`public/packages/Tutorials/bdd/features/howto/14-files/`.

The files are in Browse > Files: Demo holds the sample files that come with Datagrok, My files is the account's own
folder on the server. Opening a demo file as a table, a file from your computer, or pasted text is in
01-getting-started/02 (Open your data), and saving a Diff Studio model to My files is in 11-compute/06; these clips
don't repeat them. Everything a clip makes in My files is deleted when it ends.

## Plan

One idea per clip, most useful first.

1. Preview any file from Browse: one click shows a CSV, an image, a PDF, an SDF and a PDB file in place.
2. Upload a file to your files: drop a file from the computer onto My files, then preview and open it.
3. Organize your files: make a folder, drag a file into it, rename the file.
4. Share a folder: Share folder..., a colleague, the access level; the folder's share shows under Files.
5. Search files: the folder view's search, with or without subfolders.
6. A file's details and actions: the context panel for a file, and Download.
7. Write a text file: Create file... in My files, type in the editor, SAVE.

## 01 - Preview any file from Browse (1:17)

Title: Preview files in Datagrok: tables, images, PDFs, molecules, 3D

Description: One click on a file in Browse shows it in the main area without opening it. The video goes through
Browse > Files > Demo: a CSV shows as a table, a JPG as an image, a PDF page by page, an SDF file as a table of
drawn structures, and a PDB file as a 3D protein structure. To open a file for work, double-click it (see Open your
data in the Getting started playlist).

Feature: `public/packages/Tutorials/bdd/features/howto/14-files/01-preview-any-file.feature`

## 02 - Upload a file to your files (0:36)

Title: Upload a file to your Datagrok files

Description: My files (Browse > Files > My files) is your own folder on the Datagrok server. To upload a file, drag
it from your computer onto the folder; it shows up in the folder and in the Browse tree, where one click previews
it and a double-click opens it as a table. The video uploads a small CSV of solubility measurements.

Feature: `public/packages/Tutorials/bdd/features/howto/14-files/02-upload-a-file-to-your-files.feature`

## 03 - Organize your files (0:59)

Title: Make folders, move and rename files in Datagrok

Description: Your files in Browse > Files > My files can be kept in folders. Right-click the folder view and pick
Create folder... to make one, drag a file onto a folder in the Browse tree to move it there, and right-click a file
and pick Rename... to give it a new name. The video works on a small CSV of solubility measurements.

Feature: `public/packages/Tutorials/bdd/features/howto/14-files/03-organize-your-files.feature`

## 04 - Share a folder (1:04)

Title: Share a folder of your files in Datagrok

Description: Files are shared by folder, not one by one. Right-click a folder in Browse > Files > My files, pick
Share folder..., choose a colleague or a group and what they may do (View and use lets them read the files; Full
access adds Edit, Delete and Share), and click OK; Send notifications is optional (off in this video). The
shared folder then shows under Files, and Share folder... on it lists who has access.

Feature: `public/packages/Tutorials/bdd/features/howto/14-files/04-share-a-folder.feature`

## 05 - Search files (0:37)

Title: Find a file by name in Datagrok

Description: Every folder in Browse > Files has a search box at the top of its view. Typing part of a name finds
the matching files in the folder and in all the folders inside it; switching off the Search in subfolders icon next
to the box limits the search to the folder itself. The video finds the aspirin files in Files > Demo and opens one of
them with a double-click.

Feature: `public/packages/Tutorials/bdd/features/howto/14-files/05-search-files.feature`

## 06 - A file's details and actions (0:38)

Title: See a file's details and download it in Datagrok

Description: Click a file in Browse and open the context panel (the toggle at the bottom right of the window) to
see the file's details: the file share it is on, its folder, size and dates. The Actions section lists what you
can do with it, such as Open, Download, Rename... and Delete...; Download saves a copy to your computer. The video
uses cars.csv from Files > Demo.

Feature: `public/packages/Tutorials/bdd/features/howto/14-files/06-a-file-in-the-context-panel.feature`

## 07 - Write a text file (0:45)

Title: Write and save a text file in Datagrok

Description: Notes can live next to your data. In Browse > Files > My files, right-click the folder view and pick
Create file..., give the file a name, then click it: a text file opens in an editor. Type, and a SAVE button shows
at the top; it writes the text to your files on the server.

Feature: `public/packages/Tutorials/bdd/features/howto/14-files/07-write-a-text-file.feature`
