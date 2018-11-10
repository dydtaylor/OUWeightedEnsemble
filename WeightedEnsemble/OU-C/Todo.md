# TODO

* See all places in code with comment `// JUN COMMENT`

* Looping through time in the dynamics engine: This doesn't make sense to me. There should be a time loop inside dynamicsEngine() that loops from 0 (nt=0) to tau (t=dt*tau). There should be no dt loop in main(), right? There is a time-loop here up to tau AND a loop in dynamicsEngine? Should only be a loop in dynamicsEngine. The current code is quadratically too much?

* Reweighting after flux out: I don't understand why this while loop is necessary. Find total weight in one for-loop, then adjust everyone's weight in a second for-loop.. Many of the operations, and some of the variables, seem wasteful.

* Variable being declared in the middle of the execution block. simMaxHolder

* All loop counters - make sure they have descriptive names. Declare these, and all variables, at top of functions. The reason is that variable declaration takes cpu time, and many of these counters can be reused. Eg., you probably only need two replica counters, iReplica and jReplica.

* Use += notation

* Use isnan() instead of A!=A?

* Change the name of tau1 to tauInitialTransient or tauFirstQuarter

* Change name of tauM to tauMax.

* Since it is 1d OU, can you not define Z2, and not compute Z2? This is in a loop, in a loop, in a loop, so every microsecond counts here.
