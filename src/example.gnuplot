# example.gnuplut : configuration for plotting (change as needed)

reset                                   # reset
set size ratio 0.2                      # set relative size of plots
set grid xtics ytics                    # grid: enable both x and y lines
set grid lt 1 lc rgb '#cccccc' lw 1     # grid: thin gray lines
set multiplot layout 3,1 scale 1.0,1.0  # set two plots for this figure

# time domain
#set ylabel 'Sample value Block 10'               # set y-axis label
#set xlabel 'Sample #'                   # set x-axis label
#set yrange [-1:1]                       # set y plot range
#set xrange [0:40]                      # set x plot range
#plot '../data/ncoOut10_psd.dat' using 1:2 with lines lt 1 lw 2 lc rgb '#000088' notitle

# time domain
#set ylabel 'Sample value Blcok 11'               # set y-axis label
#set xlabel 'Sample #'                   # set x-axis label
#set yrange [-1:1]                       # set y plot range
#set xrange [0:40]                      # set x plot range
#plot '../data/ncoOut11_psd.dat' using 1:2 with lines lt 1 lw 2 lc rgb '#000088' notitle

# freq domain (Fourier)
#set ylabel 'Spectrum (Mag)'              # set y-axis label
#set xlabel 'Frequency bin'               # set x-axis label
#set yrange [-60:0]                    # set y plot range
#set xrange [0:120]                       # set x plot range
#plot '../data/fmdemod_psd.dat' using 1:2 with lines lt 1 lw 2 lc rgb '#008800' notitle

set ylabel 'Spectrum (Mag)'              # set y-axis label
set xlabel 'Frequency bin'               # set x-axis label
set yrange [-60:0]                    # set y plot range
set xrange [0:120]                       # set x plot range
plot '../data/fmstereo_psd.dat' using 1:2 with lines lt 1 lw 2 lc rgb '#008800' notitle

set ylabel 'Spectrum (Mag)'              # set y-axis label
set xlabel 'Frequency bin'               # set x-axis label
set yrange [-60:0]                    # set y plot range
set xrange [0:120]                       # set x plot range
plot '../data/fmpilot_psd.dat' using 1:2 with lines lt 1 lw 2 lc rgb '#008800' notitle

# freq domain (Fourier)
#set ylabel 'Spectrum (Mag) (Left)'              # set y-axis label
#set xlabel 'Frequency bin'               # set x-axis label
#set yrange [-145:-15]                      # set y plot range
#set xrange [0:25]                       # set x plot range
#plot '../data/stereoleft_psd.dat' using 1:2 with lines lt 1 lw 2 lc rgb '#008800' notitle

# freq domain (PSD)
set ylabel 'Spectrum (dB/Hz) (Right)'            # set y-axis label
set xlabel 'Frequency (KHz)'             # set x-axis label
set yrange [-145:-15]                       # set y plot range
set xrange [0:25]                       # set x plot range
# add your own .dat file for PSD as part of the take-home
plot '../data/stereoright_psd.dat' using 1:2 with lines lt 1 lw 3 lc rgb '#880000' notitle

unset multiplot
