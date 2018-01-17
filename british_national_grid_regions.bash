#! /bin/bash

#british_national_grid_regions.bash - get a list of squares needed from BNG
#
# Usage: extract-os-terrain.bash llx lly urx ury
# Coordinates in british national grid EPSG:27700
# Outputs a list of two-letter two-digit square names in British national grid 
#   that is needed to cover the area of a rectangle described by input 
#   coordinates to standard output.
#
# Example of usage:
#   To get a list of filenames of ascii-files e.g. SV23.asc
#
#   for reg in $(british_national_grid_regions.bash $1 $2 $3 $4)
#   do
#	echo $reg.asc
#   done
#
# Johan Arvelius, SMHI, 2017-01-17

# Return a two-letter code of a 100km square
# parameters x y
function mainregion()
{
    letters=(A B C D E F G H J K L M N O P Q R S T U V W X Y Z)
    let subx=$1
    let suby=$2
    if [ "$1" -gt 500000 ]
    then
	__firstletter=T
	if [ "$2" -gt 500000 ]
	then
	    echo "You're in the middle of the North Sea"
	    return 2
	fi
    elif [ "$2" -gt 1000000 ]
    then
	__firstletter=H
    elif [ "$2" -gt 500000 ]
    then
	__firstletter=N
    else
	__firstletter=S
    fi
    
    let __row=suby/100000
    let __row=4-__row
    let __number=__row*5
    let __number+=subx/100000 
    echo $__firstletter${letters[$__number]}
}

# Return a two-digit number of a 10km tile 
# parameters x y
function subregion()
{
    let subx=$1%100000
    let myx=$subx/10000
    let suby=$2%100000
    let myy=$suby/10000
    echo $myx$myy
}

# Return a two-letter two-digit name of a tile
# parameters x y
function region()
{
    echo $(mainregion $1 $2)$(subregion $1 $2)
}

# Return a list of tilenames needed to cover a rectange
# parameters  llx lly urx ury
function regionlist()
{
    let __llcol=$1/10000
    let __llrow=$2/10000
    let __urcol=$3/10000
    let __urrow=$4/10000
    let __nrows=$__urrow-$__llrow
    let __ncols=$__urcol-$__llcol
    for (( ii=0; ii<=$__ncols; ii++ ))
    do
	let xcoord=$1+$ii*10000
	for (( jj=0; jj<=$__nrows; jj++ ))
	do
	    let ycoord=$2+$jj*10000
	    echo $(region $xcoord $ycoord)
	done
    done
}

#mainregion=$(region $1 $2)
#echo $mainregion
#subregion=$(subregion $1 $2)
#echo $subregion
#region=$(region $1 $2)
#echo $region

#========#
## MAIN ##
#========#
regionlist $1 $2 $3 $4
