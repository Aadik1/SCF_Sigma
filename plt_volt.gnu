


set title 'Current vs Voltage, Hub = 0.1 eV'
set xlabel 'Voltge (eV)'
set ylabel 'Current (arb. units)'
set key inside bottom right
set xrange[-7:3]

#plot 'Volt_Current_0.dat' u 1:2 w l title 'No Interactions', 'Volt_Current_1_sf.dat' u 1:2 w l title 'Hartree-Fock, Hub = 1 eV'

#plot 'Volt_Current_0.dat' u 1:2 w l title 'No Interactions', 'Volt_Current_1_sf.dat' u 1:2 w l title 'Hartree-Fock, Hub = 1 eV', 'Volt_Current_2_sf.dat' u 1:2 w l title 'Second Born Approximation, Hub = 0.1 eV'

plot 'Volt_Current_2.dat' u 1:2 w l title 'Non Self-Consistent', 'Volt_Current_2_sf._H0p1dat' u 1:2 w l title 'Self-Consistent'