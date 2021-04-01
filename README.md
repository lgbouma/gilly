# gilly

This repository contains tools for gyrochronological age-dating.  The question
that spurred the development of these tools was: are there any obvious
correlations between gyrochronological age and planetary properties?

The planet samples explored included the CKS one. The rotation periods were
from Mazeh+15 (a supplement to the McQuillan ones).

Gyro-ages explored included the Angus+2019 model, and the Spaza & Landafame
2020 model.

The conclusions were "yeah, it's like 2-sigma"; these were published in
Sandoval+20 (http://arxiv.org/abs/2012.09239).  Nonetheless, these gyro
implementations are reusable.  Say you have some Gaia Bp-Rp array, and periods,
and you want the ages. Then:

```
from gilly.gyrochronology import (
  MamajekHillenbrand08_gyro,
  SpadaLanzafame20_gyro,
  Angus19_gyro
)

from cdips.utils.mamajek import get_interp_BmV_from_BpmRp

# BpmRp is the existing array of dereddened Gaia colors
BmV = get_interp_BmV_from_BpmRp(BpmRp)

# Prot is the array of periods

mh08_age = MamajekHillenbrand08_gyro(BmV, Prot)
a19_age = Angus19_gyro(BpmRp, Prot)

# can set plot as true to write a sanity check
sl20_age = SpadaLanzafame20_gyro(BmV=BmV, Prot=Prot, makeplot=False)
```
