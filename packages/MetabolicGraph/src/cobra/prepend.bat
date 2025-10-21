@echo off
set "file=sampler.js"
set "tempfile=%file%.tmp"

(
    echo export  
    type "%file%"
) > "%tempfile%"

move /Y "%tempfile%" "%file%"
