
Content
=================

<!--ts-->
   * [Timelines](#Timelines)
      * [week1](#Week1)
      * [week2](#Week2)
      * [week3](#Week3)
      * [week4](#Week4)
      * [week5_6_7](#Week5_6_7)
      * [week8](#Week8)
      * [week9](#Week9)
      * [week10_11](#Week10_11)
      * [week12](#Week12)
<!--te-->


# Project description


# Timelines

## Week1:

Read paper and learnt associated knowledge

* [Linear Transformations and the k-Means Clustering Algorithm Applications to Clustering Curves](https://github.com/sakuramomo1005/Functional-data-analysis-draft/blob/master/Wrapup0418/papers/Linear%20Transformations%20and%20the%20k-Means%20Clustering%20Algorithm%20Applications%20to%20Clustering%20Curves.pdf)
* [Stratified Psychiatry via Convexity-Based Clustering with Applications Towards Moderator Analysis](https://github.com/sakuramomo1005/Functional-data-analysis-draft/blob/master/Wrapup0418/papers/Stratified%20Psychiatry%20via%20Convexity-Based%20Clustering%20with%20Applications%20Towards%20Moderator%20Analysis.pdf)
* [pre presentation](https://github.com/sakuramomo1005/Functional-data-analysis-draft/blob/master/Wrapup0418/papers/talkFDNY.pdf)

And asked questions 
* [File](https://github.com/sakuramomo1005/Functional-data-analysis-draft/blob/master/Wrapup0418/results/some%20understandings%20and%20questions-Kate-2019-01-24.pdf)
* [Codes for this file 1](https://github.com/sakuramomo1005/Functional-data-analysis-draft/blob/master/Wrapup0418/codes/some%20understandings%20and%20questions-Kate-2019-01-24.Rmd)
* [Codes for this file 2](https://github.com/sakuramomo1005/Functional-data-analysis-draft/blob/master/Wrapup0418/codes/simulation%20and%20draw%20figure4%200124.R)


## Week2:

Aim: 

* Try to analyze the real data. Choose 2 or more baseline covariates and do k-means and see whether we could get a good VI.

The data used: 

* *longFormat_score17*

The codes: 

* [codes files](https://github.com/sakuramomo1005/Functional-data-analysis-draft/tree/master/Wrapup0418/codes/VI)

The results:

* [results files](https://github.com/sakuramomo1005/Functional-data-analysis-draft/tree/master/Wrapup0418/results/VI)

By simply choose variables without transforamtinos, the results are not very satified. 

The summary from Dr. Tarpey: [Rotation project](https://github.com/sakuramomo1005/Functional-data-analysis-draft/blob/master/Wrapup0418/results/KateRotationProject2.pdf)

## Week3:

Since the previous method did not work well. The idea got changed to the paper: [Stratified Psychiatry via Convexity-Based Clustering with Applications Towards Moderator Analysis](https://github.com/sakuramomo1005/Functional-data-analysis-draft/blob/master/Wrapup0418/papers/Stratified%20Psychiatry%20via%20Convexity-Based%20Clustering%20with%20Applications%20Towards%20Moderator%20Analysis.pdf)

Asked questions and found new focuses: 

* [summary files](https://github.com/sakuramomo1005/Functional-data-analysis-draft/tree/master/Wrapup0418/results/Purity)

## Week4:

Run Newton Raphson method to find the max purity value. 

The hcaf dataset was used. Generated a linear transformation: dat$AX = dat$BaselineCGI + dat$age

[results-week4](https://github.com/sakuramomo1005/Functional-data-analysis-draft/blob/master/Wrapup0418/results/results.pdf)

However, the results were not good at that time. There were some bugs. 

## Week5, 6, 7:

Figured out the way to calculate the purity. And drew the $\lambda$ v.s. $w$ plot

Fixed the questions on the purity 
[](https://github.com/sakuramomo1005/Functional-data-analysis-draft/blob/master/Draft/Week6/purity%20confusion%20(2).ipynb)

Began to do the simulation. 

## Week8:

Tried to use Monte Carlo method to simulate data and calculate purity. 

## Week9:

Fixed the puirty simulation method. Instead of calculating through the pdfs, generate X from standard MVN and move and shift the distributions to f1 and f2 and then calculate the purity. 

## Week10_11

* Improved previous plots by adding gamma1 and gamma2 directions. 

* Draw the similar plots with real dataset

* Try the same process with consideration of intercept: beta_0, beta_1, beta_2, and gamma_0, gamma_1, and gamma_2

* Try cvxcluster function with covariates

  * [results week 10](https://github.com/sakuramomo1005/Functional-data-analysis-draft/blob/master/Wrapup0418/results/update0331.pdf)
  * [results week 11](https://github.com/sakuramomo1005/Functional-data-analysis-draft/blob/master/Wrapup0418/results/result_0410.pdf)

## Week12

Similar simulation. But try data with 3 baseline covariates. [results week 12](https://github.com/sakuramomo1005/Functional-data-analysis-draft/blob/master/Wrapup0418/results/result_0419.pdf)

