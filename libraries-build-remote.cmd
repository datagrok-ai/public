cd js-api
cmd /c npm install
cmd /c npm run build
cd ..
forfiles /P libraries /C "cmd /c cd @file & npm install"
forfiles /P libraries /C "cmd /c cd @file & npm run build"
