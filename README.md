Regularized tau estimator
==================================

This repository contains all the code to reproduce the results of the paper
[*A new robust and efficient estimator for ill-conditioned linear inverse problems with outliers*](http://infoscience.epfl.ch/record/203430).


Abstract
--------

Solving a linear inverse problem may include difficulties such as the presence of outliers and a mixing matrix with a large condition number. In such cases a regularized robust estimator is needed. We propose a new tau type regularized robust estimator that is simultaneously highly robust against outliers, highly efficient in the presence of purely Gaussian noise, and also stable when the mixing matrix has a large condition number. We also propose an algorithm to compute the estimates, based on a regularized iterative reweighted least squares algorithm. A basic and a fast version of the algorithm are given. Finally, we test the performance of the proposed approach using numerical experiments and compare it with other estimators. Our estimator provides superior robustness, even up to 40% of outliers, while at the same time performing quite close to the optimal maximum likelihood estimator in the outlier-free case. 


Authors
-------

Marta Martinez-Camara and Martin Vetterli are with the
Laboratory for Audiovisual Communications ([LCAV](http://lcav.epfl.ch)) at 
[EPFL](http://www.epfl.ch).

<img src="http://lcav.epfl.ch/files/content/sites/lcav/files/images/Home/LCAV_anim_200.gif">

Michael Muma and Abdelhak M. Zoubir are with the
Signal Processing Group ([SPG](http://www.spg.tu-darmstadt.de/spg)) at 
[TUDarmstadt](http://www.tu-darmstadt.de).

<img src="http://www.spg.tu-darmstadt.de/media/spg/pics/spg_logo_screen_medium~1_182x0.png">


#### Contact

[Marta Martinez-Camara](mailto:marta[dot]martinez-camara[at]epfl[dot]ch) <br>
EPFL-IC-LCAV <br>
BC Building <br>
Station 14 <br>
1015 Lausanne


License
-------

Copyright (c) 2015, Marta Martinez-Camara, Michael Muma, Abdelhak M. Zoubir, Martin Vetterli

This code is free to reuse for non-commercial purpose such as academic or
educational. For any other use, please contact the authors.

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br /><span xmlns:dct="http://purl.org/dc/terms/" property="dct:title">Regularized tau estimator</span> by <a xmlns:cc="http://creativecommons.org/ns#" href="http://lcav.epfl.ch" property="cc:attributionName" rel="cc:attributionURL">Marta Martinez-Camara, Michael Muma, Abdelhak M. Zoubir, Martin Vetterli</a> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.<br />Based on a work at <a xmlns:dct="http://purl.org/dc/terms/" href="https://github.com/LCAV/RegularizedTauEstimator" rel="dct:source">https://github.com/LCAV/AcousticRakeReceiver</a>.