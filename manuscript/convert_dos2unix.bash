#!/bin/bash

perl -pe 's/\r\n|\n|\r/\n/g' $1 > $2
