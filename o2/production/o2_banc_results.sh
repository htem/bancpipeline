#!/bin/bash

# Time and start
start=`date +%s`

# Define source and destination paths
A='/n/data1/hms/neurobio/wilson/banc'
B='/n/files/Neurobio/wilsonlab/banc'

# Check if directories exist
if [ ! -d "$A" ] || [ ! -d "$B" ]; then
    echo "Error: One or both of the specified directories do not exist."
    echo "A: $A"
    echo "B: $B"
    exit 1
fi

# Hardcoded file extensions
EXTENSIONS=("png")

# Function to remove empty directories
remove_empty_dirs() {
    local dir=$1
    find "$dir" -type d -empty -delete
}

# Function to move files to 'old' directory if they don't exist in source
move_to_old_if_not_in_source() {
    local src=$1
    local dest=$2
    local old_dir="$B/matching/mirror/old"
    
    if [ -d "$dest" ]; then
        for ext in "${EXTENSIONS[@]}"; do
            find "$dest" -type f -name "*.$ext" | while read file; do
                # Skip files in directories with 'old' in the name
                if [[ "$file" == *"/old/"* || "$file" == *"/old" ]]; then
                    continue
                fi
                rel_path=${file#$dest}
                rel_path=${rel_path#/}  # Remove leading slash if present
                if [ ! -f "$src/$rel_path" ]; then
                    mkdir -p "$old_dir/$(dirname "$rel_path")"
                    mv "$file" "$old_dir/$rel_path"
                    echo "Moved $file to $old_dir/$rel_path"
                fi
            done
        done
        
        # Remove empty directories in dest and the old directory
        remove_empty_dirs "$dest"
        remove_empty_dirs "$old_dir"
        echo "Removed empty directories in $dest and $old_dir"
    else
        echo "Destination directory $dest does not exist."
    fi
}

# Function to sync files
# Function to sync files
sync_files() {
    local src=$1
    local dest=$2
    
    # Create an array of include patterns
    include_patterns=()
    for ext in "${EXTENSIONS[@]}"; do
        include_patterns+=("--include=*.$ext")
    done

    # Change permissions on source files
    find "$src" -type f \( -name "*.png" -o -name "*.csv" \) -exec chmod 644 {} +
    find "$src" -type d -exec chmod 755 {} +

    rsync -rlptv --recursive \
        --include='*/' \
        "${include_patterns[@]}" \
        --exclude='*' \
        --prune-empty-dirs \
        "$src" "$dest"

    echo "Synced files from $src to $dest"
}

# Move now outdates files on n/files to the 'old' folder, to consider for deletion
move_to_old_if_not_in_source "$A/meta/" "$B/meta/"
move_to_old_if_not_in_source "$A/matching/mirror/images/" "$B/matching/mirror/images/"
move_to_old_if_not_in_source "$A/matching/fafb/images/" "$B/matching/fafb/images/"

# Sync new files from o2
sync_files "$A/meta/" "$B/meta/"
sync_files "$A/matching/mirror/images/" "$B/matching/mirror/images/"
sync_files "$A/matching/fafb/images/" "$B/matching/fafb/images/"

# Update O2 with our manual matching results
move_to_old_if_not_in_source "$B/matching/mirror/correct/" "$A/matching/mirror/correct/"
sync_files "$B/matching/mirror/correct/" "$A/matching/mirror/correct/"
sync_files "$B/matching/fafb/correct/" "$A/matching/fafb/correct/"

# Remove any empty directories
remove_empty_dirs "$B/meta/"
remove_empty_dirs "$B/matching/mirror/images/"
remove_empty_dirs "$A/matching/mirror/correct/"
remove_empty_dirs "$B/matching/fafb/images/"
remove_empty_dirs "$A/matching/fafb/correct/"

# Anounce
echo "File transfer and movement completed."

end=`date +%s`
runtime=$((end-start))
echo "script completed in: "
echo $runtime

