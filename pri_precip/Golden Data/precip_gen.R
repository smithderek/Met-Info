---
title: "Precip Golden Data"
author: "Derek Smith"
date: "Tuesday, December 29, 2015"
output: html_document
---

```{r,echo=FALSE}
seed<-1234
set.seed(seed)
baseline<-42


```

```{r,echo=FALSE}
freq<-function(lam){if(lam=="low"){
 precip_freq<<-rpois(n=size, lambda=0.001) 
  
}

else if(lam=="medium"){
 precip_freq<<-rpois(n=size, lambda=0.005)   
  
}

else if(lam=="high"){
precip_freq<<-rpois(n=size, lambda=0.01) 
  
}

else(return("Invalid Precip Rate"))}
```

```{r,echo=FALSE}
freq(lam="low")
```

#Overview of precipitation golden data

The golden data set was created with the following inputs

* A seed `r seed`.
* An initial bucket starting point of `r baseline` mm.
* A `r lam` intensity rain event.
* Gauge noise of; `r gaugeNoise1`