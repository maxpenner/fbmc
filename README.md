# Filter Bank Multicarrier
This repository contains Filter Bank Multicarrier (FBMC) simulations for modulator, demodulator, frequency synchronization, doubly-selective channels etc.

# Scripts
- **sim_single_packet_A_nosync.m**: Generate a single frame, pass through AWGN or doubly-selective channel and measure BEP, BLEP etc.
- **sim_single_packet_B_fsync.m**: Generate a single frame, pass through AWGN or doubly-selective channel, synchronize in frequency domain and measure BEP, BLEP etc.
- **sim_model_amb_plot.m**: Plot the ambiguity function:
![ambfunc](https://user-images.githubusercontent.com/20499620/85293778-1fd35f00-b49e-11ea-8b5b-3022c3b3009d.jpg)
