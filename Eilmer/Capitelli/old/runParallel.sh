#!/bin/bash

#for dir in ~/asdf/*; do (parallel ./run-calculation-in-stages.lua :::); done

#find ~/Documents/asdf -type d | parallel "./run-calculation-in-stages.lua"

#parallel --progress ::: dir_*/run-calculation-in-stages.lua

#for dir in */; 
#do
#	echo "$dir"
#	(cd "$dir" touch asdf.txt)
#	cd "$dir"
#	parallel --jobs 3 --bar ::: ./run-calculation-in-stages.lua &
#
#	cd ..
#
#done

ls -d */ | parallel --bar --tmux --jobs 4 'cd {} && ./run-calculation-in-stages.lua'