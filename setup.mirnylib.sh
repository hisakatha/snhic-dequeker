#!/usr/bin/env bash
set -ue
#git submodule add --branch hg-python3 https://github.com/hisakatha/mirnylib-legacy
(cd mirnylib-legacy && python3 install_linux.py)
#ln -sf mirnylib-legacy/mirnylib

