# sr_ft Studentized Range - FORTRAN


1. Wrap the fortran from the _Royal Statistical Society_: (inside the statlib directory)
```
(scipydev) samuels-mbp:sr_ft swallan$ f2py -c  -m statlib statlib/cmustatlib.f
```

2. Wrap fortran from copenhaver: (missing IBM functions, will not compile)
```
samuels-mbp:sr_ft swallan$ f2py -c  -m bayreuth bayreuth/qprob.f
```


Disregard other files. 
