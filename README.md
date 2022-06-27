# SW-implementation
Python implementation of the Smith-Watermann algorithm for the course *Algorithms for bioinformatics*, held by prof. Blanzieri and prof. Tebaldi.

## Usage
The script will ask the user to input two sequences. In case of no specification of the other parameters, the default parameters will be used.

## Output
The basic output given two sequences (e.g., AAACGCT and AATCCG), will look like this:

RESULTS:
Score:  8 Match score: 3 Gap penalty: gap
SEQUENCE:  AA+C-C

Score:  8 Match score: 3 Gap penalty: gap
SEQUENCE:  AA+C+G

SEQUENCE:  AA++CG

Where *score* is the score of the alignment.
If two (or more) sequences are not separated by the *score* line (in this example, Score:  8 Match score: 3 Gap penalty: gap), it means that the sequences were reconstructued during the backtrack starting from the same initial cell.

By running the script with the -l option, a longer output containing other statistics will be printed. 

## Options
The options for the script will be displayed by running the script with -h, or can be found in the *help.txt* file.
