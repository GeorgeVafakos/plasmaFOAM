# Example Case: Varying E/N and varying gas composition

This example uses a mixture of He and artificial dry air (80% N<sub>2</sub> and 20% O<sub>2</sub>), to produce the plasma reaction rates, and the trasnport and diffusion coefficients. In this example a fixed value was given tomole fractions of the species. The reduced electric field varies from 0.1 to 100 using 30 linearly spaced points, and the mole fractions of the He and artificial dry air vary from 0 to 1.0 using 10 linearly spaced points (the air always keeps a constant proportionality of 80% N<sub>2</sub> and 20% O<sub>2</sub>). The mole fractions always sum to 1.


**Mixture:** He, N<sub>2</sub>, O<sub>2</sub> \
**Chemistry Model:** Custom model selected by the user, according to their reaction ID (./LXCatData/availableReactions.txt) \
**Reduced Electric Field:** E/N = 0.1-100 Td (linearly spaced) \
**Gas Temperature:** T = 300 K \
**Mole Fractions:** He = 0-1.0, N<sub>2</sub> = 0.8, O<sub>2</sub> = 0.2
