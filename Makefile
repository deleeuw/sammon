
SSRC = smacofIsotone.c smacofMPInverseV.c smacofSort.c smacofSSSammonEngine.c smacofSSSammonMajorize.c smacofSSSammonMonotone.c

%.o: %.c smacofSSSammon.h
	clang -c $@

shlib: smacofSSSammon.h $(SSRC)
	R CMD SHLIB -o smacofSSSammon.so $(SSRC)

clean:
	rm -f *.o

pristine:
	rm -f *.o *.so

