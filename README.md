## FISHmodel
Code used for Ietswaart et al, Cell Systems (2017) in press

Additions will be made soon.

# Stochastic cellular mRNA production and degradation simulations

Open `FCA_55.cpp` in text editor and search in `int main()`: 
choose whether you want to simulate

ON/OFF production with burst size ~ volume (scenario 1)
ON/OFF production with burst frequency ~ cell volume (scenario 2)
Poisson production (scenario 3)

by commenting out scenario 1, 2 or 3.

To run simulations: 

`chmod 755 sFCA_55.sh`

`g++ -O3 -o FCA_55 FCA_55.cpp`
 
`./sFCA_55.sh` 

