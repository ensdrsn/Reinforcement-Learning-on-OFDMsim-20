# OFDMsim-Reinforcement-Learning
Parameter optimization on OFDM simulation using Reinforcement Learning.

The purpose of this work is to use Reinforcement Learning algorithms on a communication system to 
improve transmission performance. First OFDM, the designated modulation technique 
of 5G NR, is simulated in a MATLAB environment. The created simulation works 
under various numerology configurations such as the different number of subcarriers, 
cyclic prefix (CP) lengths, modulation orders; also there are different environmental
conditions with different delay spread values. Based on this, different RL frameworks 
utilized to learn the efficient numerology configurations (actions of RL agent) under different 
environmental conditions (states of RL agent). Bit error and data rate are used as the 
performance assessment indicators of OFDM simulations.
