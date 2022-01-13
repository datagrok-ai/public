forfiles /S /M package.json /C "cmd /c npm install"
forfiles /S /M package.json /C "cmd /c npm run link-all"