  set term post eps color  enhanced dl 2 lw 2
 set encoding iso_8859_1;
 set style line 26 lt 1  lw 1  pt 7  ps 1. lc rgb    "cyan"


 set output 'L1405jl.pdf'
 set multiplot ; #  {{{

#{{{ 
 set lmargin at screen 0.2
 set rmargin at screen 0.8
 set bmargin at screen 0.95
 set tmargin at screen 0.5
 set xtics font  'Helvetica,20' offset 0,-0.5;set ytics font 'Helvetica,20' offset -0.5,0;
set ytics 4 format  '%2.0F';set mytics 2
 set xtics 0.08 format  ''; set mxtics 4; 
# set terminal png transparent nocrop enhanced font arial 8 size 420,320 
# set output 'pm3d.6.png'
set key top
 set pm3d map at s ftriangles interpolate 10,1
	 set palette defined (0 'white', 33 'green', 66 'blue',100 'red')
#	 set palette defined (0 'white', 33 'green', 66 'blue',100 'red')
#	 set palette rgbformulae 23,28,3
	 set cbtics  1 font  'Helvetica,20'
	 set cbtics 4 format  '%2.0F'
	 set xrange [1.2:1.6       ]
	 set yrange [-200:200          ]
	 set cbrange [-9:5    ]
 plot 'outputjl.txt' index 0 using 1:(($2)):($3) notitle with image
 set ylabel 'Im(z) (MeV)' font 'Helvetica,20' offset -2.9, 0.5
 unset grid 
plot 'outputjl.txt' index 0 using 1:(1000) notitle with line
 unset xlabel
 unset ylabel
 unset label
 unset ori  
#}}} 


