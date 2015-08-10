#!/bin/bash

#test LDU factorization on a large matrix
N=50; time for((k=1;k<=N;++k)) do for((i=1;i<=N;++i)) do echo -n $RANDOM " "; done; echo ; done | ./build/LDU_float
N=50; time for((k=1;k<=N;++k)) do for((i=1;i<=N;++i)) do echo -n $RANDOM " "; done; echo ; done | ./build/LDU

N=50; time for((k=1;k<=N;++k)) do for((i=1;i<=N;++i)) do echo -n \($RANDOM,$RANDOM\) " "; done; echo ; done | ./build/LDU_float
N=35; time for((k=1;k<=N;++k)) do for((i=1;i<=N;++i)) do echo -n \($RANDOM,$RANDOM\) " "; done; echo ; done | ./build/LDU
