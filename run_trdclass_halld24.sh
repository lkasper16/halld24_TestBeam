#!/bin/bash

FILEPATH=${1}
FILENAME=${2}

echo "=========>  Process FILE=$FILEPATH <=========="

root -l <<EOC
.L trdclass_halld24.C+g
trdclass_halld24 t("${FILEPATH}","${FILENAME}")
t.Loop()
EOC
