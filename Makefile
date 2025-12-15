SSRC = smacofIsotone.c smacofMPInverseV.c smacofSort.c smacofSSSammonEngine.c smacofSSSammonMajorize.c smacofSSSammonMonotone.c smacofSSElasticEngine.c smacofSSElasticMajorize.c smacofSSElasticMonotone.c

%.o: %.c smacofSSSamelas.h
	clang -c $@

shlib: smacofSSSamelas.h $(SSRC)
	R CMD SHLIB -o smacofSSSamelas.so $(SSRC)

clean:
	rm -f *.o

pristine:
	rm -f *.o *.so

