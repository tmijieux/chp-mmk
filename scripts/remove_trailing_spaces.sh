#!/bin/bash

find . \( -name '*.h' -o -name '*.c' -o -name '*.f90' \) -print0 | xargs -i -r -0 sed -r -i 's/\s*$//' {}
