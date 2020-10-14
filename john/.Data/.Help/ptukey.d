.BG
.FN Tukey
.TL
Studentized range (Tukey) distribution
.DN
Cumulative Probabilities (ptukey), Quantiles (qtukey) of Studentized Range.
More generally, these functions provide cumulative probabilities
and quantiles for the maximum of m studentized range statistics.
.CS
ptukey(q, nmeans, df, nranges=1)
qtukey(p, nmeans, df, nranges=1)
.RA
.AG q
vector of quantiles
.AG p
vector of probabilities
.AG nmeans
number of means to be compared
.AG df
degrees of freedom. If necessary this is replicated to be the same length
as p or q.
.OA
.AG nranges
In most applications this equals 1.  Otherwise this is
the number of ranges for which the maximum is taken.
.RT
cumulative probability (ptukey) or quantile (qtukey)
.SE

.DT
The studentized range is obtained by dividing the range of the set
of means by the standard error, assumed the same for all means.
Tukey's multiple range test, also known as the omega test, is based
on this statistic.

Note that eg qtukey(0.95,2,10) equals qt(0.975,10)*sqrt(2)
.SH REFERENCES
 Copenhaver, Margaret DiPonzio and Holland, Burt S: 
  Computation of the Distribution of the Maximum
  Studentized Range Statistic with Application to
  Multiple Significance Testing of Simple Effects.
  Journal of Statistical Computation and Simulation
  30, pp.1-15, 1988.
.SA

.EX
  > ptukey(3.877, 3, 10)
 0.9500129
  > qtukey(0.95, 3, 10)
 3.876777

.KW studentized range, omega test, multiple comparisons, tukey
.WR
