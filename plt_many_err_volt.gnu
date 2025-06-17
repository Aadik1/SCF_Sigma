
############################
#### Iteration Range #######

min = int(ARG1) 
max = int(ARG2)


##############################
#### Plots ###################

set title 'Error for SCF at Voltages for Hub = 1 eV'
set xlabel 'Iterations' 
set ylabel 'Error'
unset key

plot for [i=min:max] sprintf("sf_%d.dat",i) using 1:2 with lines title sprintf("Iteration:",i)

