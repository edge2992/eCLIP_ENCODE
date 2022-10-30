#!/bin/bash
cat counts.txt | awk '
{
  if($7 != 0){
    print $0
  }
}'
