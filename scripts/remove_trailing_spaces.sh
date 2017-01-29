#!/bin/bash

find . \( -name '*.h' -o -name '*.c' -o -name '*.f90' -o -name '*.[ch]pp' \) -print0 | xargs -i -r -0 sed -r -i 's/\s*$//' {}
