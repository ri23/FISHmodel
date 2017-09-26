# FISHmodel 

Code used for publication Ietswaart et al, Cell Systems 4: 622–635 (2017)
http://www.cell.com/cell-systems/fulltext/S2405-4712(17)30189-8


## Stochastic simulations of cellular *FLC* mRNA production and degradation

Open `FCA_55.cpp` in text editor and search in `int main()`: 
choose whether you want to simulate

ON/OFF production with burst size ~ volume (scenario 1)

ON/OFF production with burst frequency ~ cell volume (scenario 2)

Poisson production (scenario 3)

by uncommenting the chosen scenario and commenting out the other two.

Save these edits.

To run simulations: 

`chmod 755 sFCA_55.sh`

`g++ -O3 -o FCA_55 FCA_55.cpp`
 
`./sFCA_55.sh` 

## Stochastic simulations of *FLC* transcription for estimation of rate for mRNA release from locus

To run simulations: 

`chmod 755 sFCA_52.sh`

`g++ -O3 -o FCA_52 FCA_52.cpp`
 
`./sFCA_52.sh` 

## Cell area estimation of Z projection and cellular mRNA counting

This python script requires the JIC-CSB image analysis toolkit that can be installed by following the instructions on: 
https://github.com/JIC-CSB/FISHcount

To run area segmentation and mRNA counting: 

`python scripts/count_and_annotate.py -n 2 -r 1 -t 0.6 -s 0.8 input_data_folder/z_stack_data_file output_data_folder`

## Cell volume estimation using segmentation of Z stack

This python script requires the installation of JIC-CSB image analysis toolkit as detailed above.

To run 3D segmentation: 

`python scripts/3Dseg.py input_data_folder/z_stack_data_file --cell-level-threshold 0.6  --verbose 1`

`3Dseg.py` and `count_and_annotate.py` have been developed by Matthew Hartley and Tjelvar Olsson, affiliated with the John Innes Centre (JIC-CSB). The code is made available here with their permission.

