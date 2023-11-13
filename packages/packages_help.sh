#!/usr/bin/env bash

current_dir=$PWD

targetDir=''
if [[ $current_dir == *"/packages/"* ]]; then
    # shellcheck disable=SC2001
    targetDir="$(sed 's#.*/packages/##' <<<"$current_dir")"
    echo "Script will check links only for current package files: $targetDir"
elif [[ $PWD != *"/packages" ]]; then
    echo 'You should run script from packages directory'
    exit 1
fi

# shellcheck disable=SC2001
packages_dir="$(sed 's#packages/.*#packages#' <<<"$current_dir")"
#echo "Packages dir: $packages_dir"

targetFiles=$(grep --exclude-dir node_modules --exclude '*.js.map' --exclude-dir dist -oEr "(helpUrl|wiki) ?(:|=) ?(string = )?('|\`)(.+help/[^']+)('|\`)" . | sed -Ee "s#(.*): *.*(helpUrl|wiki) ?(:|=) ?(string = )?('|\`)(.+help\/[^']+)('|\`)#\1 \6#g")

IFS=$'\n'
for f in $(echo -e "$targetFiles"); do
    file=$(awk -F' ' '{print $1}' <<<"$f")
    link=$(awk -F' ' '{print $2}' <<<"$f")
    dir="$packages_dir/../"
    if [[ $link == "http"* ]]; then
        if ! wget -T1 -t3 --spider -q "$link" ; then
            echo "$file: helpUrl $link is unavailable"
        fi
        continue
    elif [[ $link == "\${_package.webRoot}"* ]]; then
        dir="$packages_dir/$(awk -F/ '{ print $2 }' <<<"$file")"
        # shellcheck disable=SC2001
        link=$(sed 's#${_package.webRoot}##g' <<<"$link")
    fi
    if ! [ -f "$dir/$(sed -rE -e 's,(.md)?\/?#.*,.md,g' -e 's,(.md)?$,.md,g' <<<"$link")" ]; then
        echo "$file: helpUrl $link could not be found"
    fi
done
