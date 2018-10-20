# Stochastic reachability using SReachTools

We demonstrate the features of SReachTools using various examples.

## How do I run this?

The CodeOcean capsule is setup to work right out of the box. To run any of the below scripts, go to `myMain.m` and set the appropriate flag to 1. 
Alternatively, you can see the results on the right hand side.

## What does this capsule have?

1. `SReachDynProg`: Solving the stochastic reachability of a target tube using dynamic programming
    1. `diSReachDynProg.m`: Analyze the following problems for a double integrator (**Runtime: 45 seconds**)
        1.  Stochastic viability problem
        2.  Stochastic reach-avoid problem (terminal hitting time)
        3.  Stochastic reachability of a target tube 
1. `SReachPoint`: Optimal controller synthesis for the stochastic reachability of a target tube problem. Given an initial state, we find the optimal
    controller that can drive the state of the stochastic system to stay with the target tube.
    1. `cwhSReachPoint.m` analyzes the optimal controller synthesis problem for a satellite rendezvous problem. We consider the problem of maximizing 
        the probability of the deputy satellite rendezvousing with the chief sattelite while staying within line-of-sight cone for accurate sensing
        and respecting actuation limits. The relative dynamics is given by a LTI system model (Clohessy-Wiltshire-Hill equations). 
        - Due to the absence of Gurobi, we demonstrate only `genzps-open`, `chance-open`, and `chance-affine` options, and skip `particle-open` option.
        - **Runtime: 8 minutes 51 seconds**
    1. `dubinsSReachPoint.m` analyzes the optimal controller synthesis problem for a Dubins' vehicle. We consider the problem of maximizing the 
        probability of the vehicle lying within a target tube. The vehicle has a LTV system dynamics. 
        - Due to the long time horizon, we demonstrate only `chance-open`.
        - *HSCC toolpaper example*
        - **Runtime: 29 seconds**
1. `SReachSet`: Compute the stochastic reach set for the stochastic reachability of a target tube problem. These are the set of safe initial states 
    that can be driven to stay with the target tube.
    1. `diSReachSet.m` analyzes the stochastic reachability of a target tube problem for a double integrator.
        -  Due to the low dimensionality, we can run all the options.
        - **Runtime: 6 minutes 8 seconds**
        - *HSCC toolpaper example*
    1. `cwhSReachSet.m` analyzes the stochastic reachability of a target tube problem for a satellite rendezvous problem. We consider the problem 
        identifying all the initial states from which the deputy satellite can rendezvous with the chief sattelite while staying within line-of-sight 
        cone for accurate sensing and respecting actuation limits. The relative dynamics is given by a LTI system model (Clohessy-Wiltshire-Hill 
        equations). 
        - Due to the huge time computation (about 1300 seconds) associated with `genzps-open` option, we explore only the `chance-open` option.
        Additionally, we validate the optimal controller from one of the vertices using $10^5$ Monte-Carlo simulations.
        - We obtained slight irregularities using SDPT3 as the backend solver. These irregularities were not present when using Gurobi 
        as the backend solver.
        - *HSCC toolpaper example*
        - **Runtime: 1 minute 27 seconds**
    1. `dubinsSReachSet.m` verifies a Dubins' vehicle under a known turn rate. We consider the problem of identifying the set of initial states from
        which the car can be steered to stay within a target tube under a known turn rate. 
        - Due to the long time horizon, we demonstrate only `chance-open`.
        - **Runtime: 2 minutes 15 seconds**
    