# ENSESt
ENSESt

This function uses Evolution Strategies (ES) instead of Genetic Algorithms (GA) as Evolutionary Algorithm (EA) in the NSGA-II procedure for multi-objective optimization.

The algorithm is able to find the pareto optimal front in most of the functions implemented in the file 'Examples.m', but the algorithm is unable to find the Pareto optimal front of the functions ZDT1, ZDT2, ZDT3, ZDT4 (Cases 11, 12, 13, 14 in file 'Examples.m') when the number of states is set to 30 (for ZDT1, ZDT2, ZDT3) and 10 (for ZDT4) as Prof. Deb specified in [1], but the pareto is found if the number of states is set to 2 (in those examples).

I appreciate if you are able to find the mistake in the algorithm (if there is...) and let me know to make the respective corrections.

Feel free to send me your doubts, corrections and/or suggestions. 

Thanks beforehand for downloading and reading this code, I hope it will be useful for you and other people working on multi-objective optimization.

* Bibliography:

[1] DEB, Kalyanmoy. "Multi-Objective optimization using evolutionary algorithms". John Wiley & Sons, LTD. Kanpur, India. 2004.

[2] BACK, Thomas. "Evolutionary algorithms in theory and practice". Oxford University Press. New York. 1996.

[3] BINH, Thanh; KORN, Ulrich. "MOBES: A multiobjective evolution strategy for constrained optimization problems". Third International Conference on Genetic Algorithms. Pag. 176-182. Czech Republic, 1997.

[4] BINH, Thanh. "A multiobjective evolutionary algorithm. The study cases". Technical report. Barleben, Germany. 1999.


