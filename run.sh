#!/bin/bash

if [ -d ./outputs ]; then
	rm -rf outputs
fi
mkdir outputs

perf stat -e power/energy-pkg/,power/energy-cores/,power/energy-gpu/,power/energy-ram/ ./build/fluid/fluid 2000 large.fld ./outputs/large.fld
echo "Finished execution"
