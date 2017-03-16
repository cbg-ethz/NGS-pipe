#!/bin/bash
#

grep -v 'REJECT' $1 > $3
grep -v 'REJECT' $2 > $4
