##!/bin/bash

#pushd ..
#Vendor/Binaries/Premake/macOS/premake5 --cc=clang --file=Build.lua gmake2
#popd

#!/bin/bash
set -e

# Move to project root (script is in Scripts/, so go up one level)
cd "$(dirname "$0")/.."

# Download Premake if not present
PREMAKE_PATH="Vendor/Binaries/Premake/macOS/premake5"
if [ ! -f "$PREMAKE_PATH" ]; then
    echo "Downloading Premake..."
    mkdir -p Vendor/Binaries/Premake/macOS
    curl -L https://github.com/premake/premake-core/releases/download/v5.0.0-beta2/premake-5.0.0-beta2-macosx.tar.gz \
        | tar -xz -C Vendor/Binaries/Premake/macOS
fi

# Generate Makefiles
Vendor/Binaries/Premake/macOS/premake5 --cc=clang --file=Build.lua gmake2

# Build
make config=release

BINARY=$(find Binaries -path "*/Release/App/App" | head -1)

if [ -z "$CI" ]; then
    sudo install -m 755 "$BINARY" /usr/local/bin/raytrace
    echo "Installed! To run the raytracer, try: raytrace"
else
    echo "Build successful!"
fi

