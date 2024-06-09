include makefile.defs

xcog:
	cd util; $(MAKE); cd ..
	cd curves; $(MAKE); cd ..
	cd stretchings; $(MAKE); cd ..
	cd mappings; $(MAKE); cd ..
	if [ "${TG_FLAGS}" = "-DTHE_GAR" ]; then \
          cd mappings/theodorsen_garrick; $(MAKE); cd ../..; \
        fi
	cd overlap; $(MAKE); cd ..
	cd src; $(MAKE); cd ..

clean:
	cd util; $(MAKE) clean; cd ..
	cd curves; $(MAKE) clean; cd ..
	cd stretchings; $(MAKE) clean; cd ..
	cd mappings; $(MAKE) clean; cd ..
	cd mappings/theodorsen_garrick; $(MAKE) clean; cd ../..
	cd overlap; $(MAKE) clean; cd ..
	cd src; $(MAKE) clean; cd ..



