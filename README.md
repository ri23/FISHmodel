# FISHmodel
Code used for Ietswaart et al, Cell Systems (2017) in press

Additions will be made soon.

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

## Stochastic simulations of *FLC* transcription



