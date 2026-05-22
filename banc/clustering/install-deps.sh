#!/usr/bin/env bash
# Install Python dependencies for banc-spectral-clustering.py
# Checks for each package first and only installs what's missing.

set -e

PACKAGES=(numpy pandas scipy scikit-learn umap-learn plotly pyarrow)

missing=()
for pkg in "${PACKAGES[@]}"; do
    # Map pip package names to Python import names
    case "$pkg" in
        scikit-learn) import_name="sklearn" ;;
        umap-learn)   import_name="umap" ;;
        pyarrow)      import_name="pyarrow" ;;
        *)            import_name="$pkg" ;;
    esac

    if ! python3 -c "import ${import_name}" 2>/dev/null; then
        missing+=("$pkg")
    fi
done

if [ ${#missing[@]} -eq 0 ]; then
    echo "All dependencies already installed."
else
    echo "Installing missing packages: ${missing[*]}"
    pip install "${missing[@]}"
fi
