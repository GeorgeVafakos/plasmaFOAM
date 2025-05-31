# AC-DBD Demo Case

This tutorial simulates the AC-driven sDBD asymmetric actuator, using a simple geometry, larger electrodes. The plasma actuator is placed in quiescent air and is powered by a sinusoidal high voltage $V_{p-p}=20~ kV$ and frequency $f=10~ kHz$. 

The thickness if the dielectric plate of the actuator is $3~mm$, and the electrodes' $0.5~mm$. The length of the high-voltage and the gounded electrode is $5$ and $20~mm$, respectively. There is no inter-electrode gap. The dielectric plate is made of Kapton $(\varepsilon_r=3.0)$ and the grounded electrode is encapsulated into dielectric material with dielectric constant $\varepsilon_r=3.7$.

In this demo case. 5 discharge events have been inmposed during the first half of the positive negative-cycle (between $T/2$ and $3T/4$). The distance between the discharge events is $6.25~\mu s$. The charge density of each discharge is set to $\rho_c=0.01~C/m$ 

To run the three stage mdoel, simply execute the `AllCasesRun` script. 

At the ionic motion stage (2<sup>nd</sup> stage), the python script `plotProbesPy.py` plots the charge density at specific locations, as defined in the `probes` dictionary, at the `system` directory. To execute the `plotProbesPy.py`, navigate to the `2_ionicMotion` directory and execute:

```bash
$ python3 plotProbesPy.py
```
Note: The linux distribution must have python3 installed, along with the matplotlib library (both always preinstalled).


More information of the usage of the three-stage model and the solution procedure can be found in Vafakos *et al.* [1]

## References

1. Vafakos, G.P., Papadopoulos, P.K., Svarnas, P. A three-stage plasma model based on one-way coupling of plasma dynamics, ionic motion, and fluid flow: Application to DBD plasma actuators. *J Appl Phys 137*, 043302 , (2025). https://doi.org/10.1063/5.0242676