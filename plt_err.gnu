############################
#### File Selection #######

i = int(ARG1)
j = int(ARG2)

set title 'SCF Error per Iteration with Pulay = 0.1'
set xlabel 'Iteration'
set ylabel 'Error'
set key inside top right
#set yrange[0:1]

plot sprintf("sf_%d.dat",i)with lines title sprintf("Iteration: %d ", i), sprintf("sf_%d.dat",j)with lines title sprintf("Iteration: %d", j)