# FourRussiansRNA
Implementation of the two vectors / four russians for accelerating Nussinov RNA folding prediction algorithm. The speedup makes the algorithm go from a **O(n^3)** to a **O(n^3/log(n))** time complexity.
The algorithm used is discribed in [Venkatachalam et al. paper](https://almob.biomedcentral.com/articles/10.1186/1748-7188-9-5).

The algorithm was implemented in c++. We've also implemented Nussinov algorithm in order to test the speed-up. Note that after crossvalidating we found that the best size of q si log in base 2 of n, this is thus the default size of q chosen.

## Presentation
[Here](https://github.com/YannDubs/FourRussiansRNA/blob/master/OnlyThreeRussians.pdf) you can find an intuitive explanation of the two-vectors acceleration, as well as our results.

## Code
The code can be found [here](https://github.com/YannDubs/FourRussiansRNA/tree/master/FourRussiansRNA) simply compile and run. Note that you can specify the files containing the sequences in the main().
