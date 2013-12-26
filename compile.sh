#!/bin/bash

##
# C++ and Fortran Compiler Script
#
# This script is designed to take several input files
# as arguments and utilize the GCC C++ and Fortran
# compilers to combine them into one executable file.
#
# @author    Oliver Spryn
# @copyright Copyright (c) 2013 and Onwards, Grove City College
# @namespace global
# @package   global
# @since     1.0.0
#

##
# Configuration
# --------------------------
#

cCompiler="gcc"
fortranCompiler="f95"

cExtn=("c", "cc", "cpp", "cxx", "c++", "h", "hh", "hpp", "hxx", "h++")
fExtn=("f", "for", "f90", "f95")

cleanUp=1
objExtn=".o"
output="run.exe"

black="\033[30m"
blue="\033[34m"
cyan="\033[36m"
green="\033[32m"
magenta="\033[35m"
red="\033[31m"
reset="\033[0m"
white="\033[37m"
yellow="\033[33m"

##
# Instructions
# --------------------------
#

function instructions {
	echo -en "$yellow"
	echo -e "C++ and Fortran Compiler Script$reset"
	echo
	echo "This script is designed to take several input files as arguments and utilize the GCC C++ and Fortran compilers to combine them into one executable file."
	echo
	echo "To use this script, follow the following convention:"
	echo -e "$red    ./compile.sh [OPTIONS] [ARGUMENTS]$reset"
	echo
	echo -e "Where $red[ARGUMENTS]$reset can be a C++ file with any of the following extensions:"
	echo -e "$blue    ${cExtn[@]}$reset"
	echo "and/or a Fortran file with any of the following extensions:"
	echo -e "$blue    ${fExtn[@]}$reset"
	echo
	echo -e "All of the files in the $red[ARGUMENTS]$reset will be compiled into a single executable."
	echo
	echo -e "The following are configurable $red[OPTIONS]$reset:"
	echo -e "$red    -d$reset"
	echo "      Disable the automatic clean up the temporary object files which are genereted by this script."
	echo "      Default: Cleans up temporary object files"
	echo
	echo -e "$red    -o <file>$reset"
	echo "      Generate an executable file with a name which was assigned by this flag."
	echo "      Default: run.exe"
	echo
	echo "A complete example of a valid input string would be:"
	echo -e "$red    ./compile.sh -d -o run.exe main.cpp main.f$reset"
}

##
# Compiler
# --------------------------
#

function compile {
	for f in $@
	do
		if [[ "${f:0:1}" != "-" && "$f" != "$1" ]]
		then
		# Determine which compiler to use, based on the extension
			extn=$(echo $f | rev | cut -d "." -f1 | rev)
			compiler="none"

		# Is this a C++ source file?
			for i in ${cExtn[@]}
			do
				if [[ "${i: -1:1}" == "," ]]
				then
					i="${i:0: -1}"
				fi

				if [[ "$i" == "$extn" ]]
				then
					compiler=$cCompiler
					break
				fi
			done

		# Is this a Fortran source file?
			if [[ "$compiler" == "none" ]]
			then
				for i in ${fExtn[@]}
				do
					if [[ "${i: -1:1}" == "," ]]
					then
						i="${i:0: -1}"
					fi

					if [[ "$i" == "$extn" ]]
					then
						compiler=$fortranCompiler
						break
					fi
				done
			fi

		# Is this file neither a C++ nor a Fortran source file?
			if [[ "$compiler" == "none" ]]
			then
				echo -e "$red[FATAL ERROR] Given a file with an unknown extension: $extn$reset"
				exit 1
			fi

		# Now run the appropriate compiler
			if [[ "$compiler" == "$cCompiler" ]]
			then
				echo "Compiling C/C++ source file: $f"
				$cCompiler -c -o $f$objExtn $f
			else
				echo "Compiling Fortran source file: $f"
				$fortranCompiler -c -o $f$objExtn $f
			fi

		# Did the compiler give any errors?
			if [[ $? -ne 0 ]]
			then
				echo -e "$red[FATAL ERROR] Compiler script terminated$reset"
				clean $@
				exit 1
			fi
		fi
	done
}

function link {
	cmd=""

	for f in $@
	do
		if [[ "${f:0:1}" != "-" && "$f" != "$1" ]]
		then
			cmd=$cmd" "$f$objExtn
		fi
	done

# Run the linker to combine all of the object files into an executable
	echo "Linking files:"$cmd
	$cCompiler -o $1 $cmd -lgfortran

# Did the linker give any errors?
	if [[ $? -ne 0 ]]
	then
		echo -e "$red[FATAL ERROR] Compiler script terminated$reset"
		clean $@
		exit 1
	fi
}

function clean {
	if [[ cleanUp -eq 1 ]]
	then
		for f in $@
		do
			if [[ "${f:0:1}" != "-" && "$f" != "$1" ]]
			then
				unlink $f$objExtn
			fi
		done
	fi
}

##
# Program
# --------------------------
#

clear

# Get all of the flags to the program
minArgs=0

while getopts :o:d option
do
	case $option in
		d)
			cleanUp=0
			((minArgs++));;

		o)
			output=$OPTARG
			((minArgs+=2));;

		?)
			echo -e "Unknown flag, ignoring: $red$OPTARG$reset Supported flags: $blue-d -o <file>$reset"
			((minArgs++));;
	esac
done

# Branch the script based on the number of arguments
if [[ $# -le $minArgs ]]
then
	instructions
else
# Display the name of the executable which will be generated
	echo -e "Output compiled file to: $cyan$output$reset"

# Should the script clean up temporary object files?
	if [[ cleanUp -eq 1 ]]
	then
		echo -e "Clean up temporary object files:$cyan true$reset"
		steps=3
	else
		echo -e "Clean up temporary object files:$cyan false$reset"
		steps=2
	fi

	echo

# Run the compiliers
	echo -e "$green[Step 1 of $steps] Compiling source files$reset"
	compile $output $@

# Run the linker
	echo -e "$green[Step 2 of $steps] Linking object files$reset"
	link $output $@

# Clean up object files
	if [[ cleanUp -eq 1 ]]
	then
		echo -e "$green[Step 3 of 3] Cleaning up temporary object files$reset"
		clean $output $@
		echo "Done"
	fi

	echo -e "\n\033[30;42mCompleted$reset\n"
fi
