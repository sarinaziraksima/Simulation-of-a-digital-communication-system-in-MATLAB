# Simulation-of-a-digital-communication-system-in-MATLAB

In this code, we simulated a system in which random bits are generated, and then channel coding is applied so that it could be modulated as a signal. This signal is transferred through the transmitter filter and the channel using a predefined transfer functions. This is followed by adding noise to simulate a roughly realistic channel. After the signal goes through the receiver filter, the received data is decoded and presented along with the number of incorrect bits as errors.

Note: The program is customizable, so initial values can be changed. Additionally, to check the system with no noise, you can set the SNR to 100.

To further examine the accuracy of this system, we used 1000 or more samples to generate a diagram that shows the error percentage at each SNR.
