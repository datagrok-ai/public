#!/bin/bash

#set -ex

# Input file
input_file="${1:-input.txt}"

# Check if input file exists
if [ ! -f "$input_file" ]; then
    echo "Input file '$input_file' does not exist."
    exit 1
fi

temp_file="temp.txt"
grep -e 'Broken link' -e 'Broken anchor' -e 'linking to' -e "markdown link couldn't be resolved" "$input_file" > $temp_file || true

if [ -f "${temp_file}" ]; then
    current_line_temp_file="temp_line.txt"

    output_file="${2:-output.csv}"
    echo '"source_page","link_to","resolved_to"' > "$output_file"

    # Read input file line by line
    while IFS= read -r line; do
        # If the line contains "->"
        if [[ $line =~ "Broken link on source page path" || $line =~ "Broken anchor on source page path" ]]; then
            echo "$(sed -nE "s# *- Broken (link|anchor) on source page path = ([^ ]+):#\2#p" <<<$line)" > "$current_line_temp_file"
        elif [[ $line == *"->"* ]]; then
            # Replace "->" with the content of the closest preceding line
            link=$(sed -nE "s# +-> linking to ([^ ]+) ?#\1#p" <<<$line | sed 's/([^)]*)//g')
            resolved=$(sed -nE "s#.*\(resolved as: ([^ ]+)\)#\1#p" <<<$line)
            echo "\"$(cat $current_line_temp_file)\",\"${link}\",\"${resolved:-$link}\"" >> "$output_file"
        elif [[ $line == *"markdown link couldn't be resolved"* ]]; then
            link=$(sed -nE "s#Error: Docs markdown link couldn't be resolved: \(([^ ]+)\) in \"([^ ]+)/(help/[^ ]+)\" for version .*#\"\1\",\"\3\",\"\"#p" <<<$line)
            echo "$link" >> "$output_file"
        else
            echo "$line" >> "$output_file"
        fi
    done < "$temp_file"

    rm -f "$temp_file" >/dev/null 2>&1
    rm -f "$current_line_temp_file" >/dev/null 2>&1
else
    echo 'No broken link information found'
fi
