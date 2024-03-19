# UltrasoundRP
The scripts for neuron signal processing and visualization, and the generation and validation of arbitrary acoustic patterns.

Scripts were tested in MATLAB 2019b.
AS_APM_256.m is the script using the angular spectrum method to compute the amplitude and phase distribution of the 2D array to generate arbitrary patterns. A "V.png" was included for testing. In the end of the script, a direct integral method was also included to validate the pattern generation. 
NeuronDataProcessing.m is used to visualize, process, and analyze neuron signals received by the Lablynx system. A sample dataset is shared at https://figshare.com/articles/dataset/UltrasoundRP/25438252 for testing.

The running time of AS_APM_256.m is within one second. It can be significantly longer if more image or simulation pixels were used. 
The running time of NeuronDataProcessing.m is around 10 to 15 minutes, depending on the CPU of the testing platform and the number of threads in parallel computing.

The test platform: Intel i7-7700K@3.60 GHz, Windows 11 and MATLAB 2019b. 
