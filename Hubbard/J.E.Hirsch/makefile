soc=Hubbard.f90
mod=band
$(mod):$(soc)
	f2py -c $(soc) -m $(mod)
clean:
	rm -f $(mod).so
