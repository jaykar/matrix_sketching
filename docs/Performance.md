# Performance

## Boost vs. Armadillo vs. Intel
![Instantiation](./images/instantiation.png)
Boost is a several times slower (3-10x), while armadillo and intel are equally fast
![Rand Fill](./images/rand_fill.png)
Interestingly, Intel performed the worst, followed by boost then Armadillo. This is because Armadillo's library rand_n function is optimized more so
![Subtract](./images/subtract.png)
Just like in instantiation, Boost trailed while the other two remained similar
![Multiplication](./images/multiplication.png)
In multiplication, the boost implementation of prod() clearly lags behind Intel and Armadillo's 
![Armadillo vs. Intel Mult](./images/arm_intel.png)
With Boost's data removed, one can see that Intel MKL is vastly faster than Armadillo's
![Concatenation](./images/concat.png)
Like before, Boost is slowe (7x), with Intel edging Armadillo slightly