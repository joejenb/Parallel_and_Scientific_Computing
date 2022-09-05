## Overview 

A C++ implementation of the colliding suns problem, in which collisions result in the formation of a new body if above a certain value of mass and velocity. OpenMP is also used to varying extents in each stage to facilitate vectorisation and parallelisation.

## Step 1

Baseline implementation

## Step 2

Added in vectorisation for update of acceleration

## Step 3

Extended calculations to make use of Runge-Kutta(2) method, hence increasing the precision of calculations.

## Step 4

Added in parallelisation of calculations.

## Speed Up

![Screenshot from 2022-09-05 14-45-27](https://user-images.githubusercontent.com/30124151/188468485-edc60b94-4e0a-452d-9c15-e3f791e93507.png)
