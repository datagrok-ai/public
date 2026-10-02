@echo off
set "file=%~dp0sampler.js"
set "tempfile=%file%.tmp"

(
    echo export  
    type "%file%"
) > "%tempfile%"

move /Y "%tempfile%" "%file%"
