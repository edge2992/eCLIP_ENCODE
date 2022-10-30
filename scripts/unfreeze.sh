#!/bin/bash
#
find ./data/eCLIP -type f -name "*.gz" -exec gunzip {} \;
