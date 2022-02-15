cd js-api
rd /S /Q node_modules
cmd /c npm install
cmd /c npm link
cmd /c npm run build
cd ..
forfiles /P libraries /C "cmd /c cd @file & rd /S /Q node_modules"
forfiles /P libraries /C "cmd /c cd @file & npm install"
forfiles /P libraries /C "cmd /c cd @file & npm link"
forfiles /P libraries /C "cmd /c cd @file & npm run link-all"
forfiles /P libraries /C "cmd /c cd @file & npm run build"
